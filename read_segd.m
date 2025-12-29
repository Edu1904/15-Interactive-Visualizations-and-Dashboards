function segd = read_segd(filename, varargin)
%READ_SEGD Read SEG-D seismic data files in MATLAB.
%   segd = READ_SEGD(filename) reads a demultiplexed SEG-D file using big
%   endian byte order and returns a structure containing the raw headers,
%   decoded metadata, and trace samples as a double matrix of size
%   [samplesPerTrace x traceCount].
%
%   Name-value pairs allow you to describe the layout of your data when the
%   header information is insufficient or proprietary:
%     'TraceCount'       - Number of traces in the file/record.
%     'SamplesPerTrace'  - Samples per trace (required if TraceCount is not
%                          provided).
%     'TraceHeaderBytes' - Bytes to skip before each trace (default 0).
%     'SampleFormat'     - Storage type for samples. Supported values:
%                          'int16', 'int24', 'int32', 'single' (default
%                          'int32'). Data are returned as double.
%     'DataOffset'       - Byte offset for the first trace. When omitted,
%                          it is estimated from the general/external
%                          header counts.
%     'Preset'           - Convenience preset; supports 'sercel_wing'
%                          (32-bit IEEE, Sercel-specific headers).
%     'GeneralHeaderBlocks'    - Override count of 32-byte general headers.
%     'ExternalHeaderBlocks'   - Override count of 32-byte external headers.
%     'AdditionalHeaderBlocks' - Override count of 32-byte additional headers.
%     'ScanTypeHeaderBytes'    - Override scan type header size (bytes).
%     'ExtendedHeaderBytes'    - Override extended header size (bytes).
%     'ExternalHeaderBytes'    - Override external header size (bytes).
%
%   Example:
%     segd = read_segd('line01.sgd', ...
%         'SamplesPerTrace', 4000, ...
%         'TraceHeaderBytes', 20, ...
%         'SampleFormat', 'int32');
%
%   The function is intentionally lightweight; SEG-D variants differ, so
%   supply explicit layout information when your file deviates from the
%   standard block sizes.

arguments
    filename (1, :) char
end

arguments (Repeating)
    varargin
end

p = inputParser;
p.addParameter('TraceCount', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.addParameter('SamplesPerTrace', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.addParameter('TraceHeaderBytes', 0, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('SampleFormat', 'int32', @(x) ischar(x) || isstring(x));
p.addParameter('DataOffset', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
p.addParameter('Preset', '', @(x) ischar(x) || isstring(x));
p.addParameter('GeneralHeaderBlocks', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
p.addParameter('ExternalHeaderBlocks', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
p.addParameter('AdditionalHeaderBlocks', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
p.addParameter('ScanTypeHeaderBytes', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
p.addParameter('ExtendedHeaderBytes', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
p.addParameter('ExternalHeaderBytes', [], @(x) isempty(x) || (isscalar(x) && x >= 0));
p.parse(varargin{:});
opts = p.Results;

preset = lower(string(opts.Preset));
if preset == "sercel_wing"
    if strcmpi(opts.SampleFormat, 'int32')
        opts.SampleFormat = 'single'; % format code 8058 => 32-bit IEEE
    end
    if isempty(opts.GeneralHeaderBlocks), opts.GeneralHeaderBlocks = 3; end
    if isempty(opts.AdditionalHeaderBlocks), opts.AdditionalHeaderBlocks = 2; end
    if isempty(opts.ScanTypeHeaderBytes), opts.ScanTypeHeaderBytes = 512; end
    if isempty(opts.ExtendedHeaderBytes), opts.ExtendedHeaderBytes = 1024; end
    if isempty(opts.ExternalHeaderBytes), opts.ExternalHeaderBytes = 1024; end
end

sampleFormat = lower(string(opts.SampleFormat));
switch sampleFormat
    case "int16"
        sampleBytes = 2;
        freadType = '*int16';
    case "int32"
        sampleBytes = 4;
        freadType = '*int32';
    case {"single", "float32"}
        sampleBytes = 4;
        freadType = '*single';
        sampleFormat = "single";
    case "int24"
        sampleBytes = 3;
        freadType = 'uint8=>uint8'; % handled manually
    otherwise
        error("Unsupported SampleFormat '%s'.", sampleFormat);
end

[fid, msg] = fopen(filename, 'r', 'ieee-be');
if fid == -1
    error('Unable to open "%s": %s', filename, msg);
end
cleanup = onCleanup(@() fclose(fid));

fileInfo = dir(filename);
if isempty(fileInfo) || fileInfo.bytes == 0
    error('File "%s" is empty or not accessible.', filename);
end
fileSize = double(fileInfo.bytes);

% Read the mandatory first 32-byte general header block.
gh1 = fread(fid, 32, 'uint8=>uint8');
if numel(gh1) < 32
    error('File "%s" does not contain a complete SEG-D general header.', filename);
end
formatCode = bitor(bitshift(uint16(gh1(3)), 8), uint16(gh1(4)));
baseScanCode = gh1(23);
switch baseScanCode
    case 4
        baseScanMs = 0.25;
    case 8
        baseScanMs = 0.5;
    case 10
        baseScanMs = 1;
    case 20
        baseScanMs = 2;
    case 40
        baseScanMs = 4;
    otherwise
        baseScanMs = NaN;
end

generalHeaderBlocks = max(1, bcd2dec(gh1(27)));
externalHeaderBlocks = bcd2dec(gh1(29));
additionalHeaderBlocks = bcd2dec(gh1(30));
if ~isempty(opts.GeneralHeaderBlocks), generalHeaderBlocks = opts.GeneralHeaderBlocks; end
if ~isempty(opts.ExternalHeaderBlocks), externalHeaderBlocks = opts.ExternalHeaderBlocks; end
if ~isempty(opts.AdditionalHeaderBlocks), additionalHeaderBlocks = opts.AdditionalHeaderBlocks; end
scanTypeHeaderBytes = opts.ScanTypeHeaderBytes;
extendedHeaderBytes = opts.ExtendedHeaderBytes;
externalHeaderBytes = opts.ExternalHeaderBytes;

if isempty(opts.DataOffset)
    headerBytes = 32 * (generalHeaderBlocks + additionalHeaderBlocks);
    if ~isempty(externalHeaderBlocks)
        headerBytes = headerBytes + 32 * externalHeaderBlocks;
    end
    if ~isempty(scanTypeHeaderBytes)
        headerBytes = headerBytes + scanTypeHeaderBytes;
    end
    if ~isempty(extendedHeaderBytes)
        headerBytes = headerBytes + extendedHeaderBytes;
    end
    if ~isempty(externalHeaderBytes)
        headerBytes = headerBytes + externalHeaderBytes;
    end
    dataOffset = headerBytes;
    % Guard against malformed headers; fall back to the first header block.
    if dataOffset <= 0 || dataOffset >= fileSize
        dataOffset = 32 * generalHeaderBlocks;
    end
else
    dataOffset = opts.DataOffset;
end

dataOffset = min(max(0, round(dataOffset)), fileSize);
fseek(fid, 0, 'bof');
rawHeader = fread(fid, dataOffset, 'uint8=>uint8');
if numel(rawHeader) < dataOffset
    error('Unable to read the declared header bytes from "%s".', filename);
end
fseek(fid, dataOffset, 'bof');

availableBytes = fileSize - dataOffset;
if availableBytes <= 0
    error('No data found after the headers in "%s".', filename);
end

if isempty(opts.TraceCount) && isempty(opts.SamplesPerTrace)
    error('Provide either TraceCount or SamplesPerTrace so the layout can be derived.');
end

traceHeaderBytes = double(opts.TraceHeaderBytes);
traceCount = opts.TraceCount;
samplesPerTrace = opts.SamplesPerTrace;

if isempty(samplesPerTrace) && ~isempty(traceCount)
    samplesPerTrace = floor((availableBytes - traceCount * traceHeaderBytes) / (sampleBytes * traceCount));
elseif isempty(traceCount) && ~isempty(samplesPerTrace)
    traceCount = floor(availableBytes / (traceHeaderBytes + samplesPerTrace * sampleBytes));
end

if samplesPerTrace <= 0 || traceCount <= 0
    error('Derived layout is invalid. Check TraceCount, SamplesPerTrace, and TraceHeaderBytes.');
end

bytesPerTrace = traceHeaderBytes + samplesPerTrace * sampleBytes;
requiredBytes = bytesPerTrace * traceCount;
if requiredBytes > availableBytes
    error(['File size (%d bytes after headers) is insufficient for %d traces of %d samples. ' ...
           'Adjust TraceCount, SamplesPerTrace, TraceHeaderBytes, or DataOffset.'], ...
           availableBytes, traceCount, samplesPerTrace);
end

traces = zeros(samplesPerTrace, traceCount, 'double');
traceHeaders = cell(traceCount, 1);

for tr = 1:traceCount
    if traceHeaderBytes > 0
        traceHeaders{tr} = fread(fid, traceHeaderBytes, 'uint8=>uint8');
    else
        traceHeaders{tr} = uint8([]);
    end

    if sampleFormat == "int24"
        raw = fread(fid, samplesPerTrace * 3, freadType);
        if numel(raw) < samplesPerTrace * 3
            error('Unexpected end of file while reading samples for trace %d.', tr);
        end
        raw = reshape(raw, 3, []);
        vals = int32(raw(1, :)) * 65536 + int32(raw(2, :)) * 256 + int32(raw(3, :));
        neg = vals >= 2^23;
        vals(neg) = vals(neg) - 2^24; % sign extension for 24-bit
        traces(:, tr) = double(vals);
    else
        vals = fread(fid, samplesPerTrace, freadType);
        if numel(vals) < samplesPerTrace
            error('Unexpected end of file while reading samples for trace %d.', tr);
        end
        traces(:, tr) = double(vals);
    end
end

segd = struct( ...
    'filename', filename, ...
    'data_offset', dataOffset, ...
    'format_code', formatCode, ...
    'base_scan_interval_ms', baseScanMs, ...
    'general_header_blocks', generalHeaderBlocks, ...
    'external_header_blocks', externalHeaderBlocks, ...
    'additional_header_blocks', additionalHeaderBlocks, ...
    'scan_type_header_bytes', scanTypeHeaderBytes, ...
    'extended_header_bytes', extendedHeaderBytes, ...
    'external_header_bytes', externalHeaderBytes, ...
    'raw_header', rawHeader, ...
    'trace_header_bytes', traceHeaderBytes, ...
    'samples_per_trace', samplesPerTrace, ...
    'trace_count', traceCount, ...
    'sample_format', char(sampleFormat), ...
    'trace_headers', {traceHeaders}, ...
    'traces', traces);

end

function out = bcd2dec(bytes)
%BCD2DEC Convert packed BCD bytes to decimal.
    if isempty(bytes)
        out = [];
        return;
    end
    digits = zeros(1, numel(bytes) * 2);
    for k = 1:numel(bytes)
        hi = bitshift(bytes(k), -4);
        lo = bitand(bytes(k), 15);
        digits(2 * k - 1) = hi;
        digits(2 * k) = lo;
    end
    if any(digits > 9)
        out = [];
        return;
    end
    out = 0;
    for d = digits
        out = out * 10 + d;
    end
end
