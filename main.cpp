#include <cuda.h>

#include <stdio.h>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include <Windows.h>
#include <psapi.h>

#undef ENABLE_PIX_RUNTIME
#define ENABLE_PIX_RUNTIME 0

#if ENABLE_PIX_RUNTIME
#include <Windows.h>
#include <pix3.h>
#include <ctime>
#endif

#define K_NOOP(...) (void)sizeof(0, __VA_ARGS__)

#if ENABLE_PIX_RUNTIME
#define K_PIX_BEGIN_EVENT(colorARGB, str) PIXBeginEvent(colorARGB, str)
#define K_PIX_END_EVENT()                 PIXEndEvent()
#else
#define K_PIX_BEGIN_EVENT(colorARGB, str) K_NOOP()
#define K_PIX_END_EVENT()                 K_NOOP()
#endif

#define CUDA_SAFE_CALL(x)                                   \
	do {                                                    \
		CUresult result = x;                                \
		if(result != CUDA_SUCCESS) {                        \
			const char *msg;                                \
			cuGetErrorName(result, &msg);                   \
			fprintf(stderr,                                 \
				"\nError: " #x " failed with error '%s'\n", \
				msg);                                       \
			exit(1);                                        \
		}                                                   \
	} while(0)


template <size_t N>
inline void HumanReadableByteSizePrintf(char (&str)[N], size_t bytes)
{
	const int32_t exponent = bytes > 0 ? static_cast<int32_t>(log10(static_cast<double>(bytes)) / log10(1024)) : 0;
	snprintf(str, N, "%.3f%c%s", bytes / pow(1024.f, exponent), "BKMGTPE"[exponent], exponent > 0 ? "B" : "");
}

// Byte size format helper
// Usage:
//		Log(HumanReadableByteSizeStr(bytes));
//		Log("Byte size: %s", HumanReadableByteSizeStr(bytes).c_str());
struct HumanReadableByteSizeStr
{
	HumanReadableByteSizeStr(size_t bytes)
	{
		HumanReadableByteSizePrintf(m_text, bytes);
	}

	operator char const *() const { return m_text; }
	char const * c_str() const { return m_text; }

	private:
		char m_text[256];
};

#include <strsafe.h>

void ErrorExit(LPTSTR lpszFunction) 
{ 
	// Retrieve the system error message for the last-error code

	LPVOID lpMsgBuf;
	LPVOID lpDisplayBuf;
	DWORD dw = GetLastError(); 

	FormatMessage(
		FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		FORMAT_MESSAGE_FROM_SYSTEM |
		FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		dw,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(LPTSTR) &lpMsgBuf,
		0, NULL );

	// Display the error message and exit the process

	lpDisplayBuf = (LPVOID)LocalAlloc(LMEM_ZEROINIT, 
		(lstrlen((LPCTSTR)lpMsgBuf) + lstrlen((LPCTSTR)lpszFunction) + 40) * sizeof(TCHAR)); 
	snprintf((LPTSTR)lpDisplayBuf, 
		LocalSize(lpDisplayBuf) / sizeof(TCHAR),
		TEXT("%s failed with error %d: %s"), 
		lpszFunction, dw, (LPCTSTR)lpMsgBuf); 
	MessageBox(NULL, (LPCTSTR)lpDisplayBuf, TEXT("Error"), MB_OK); 

	LocalFree(lpMsgBuf);
	LocalFree(lpDisplayBuf);
	ExitProcess(dw); 
}

struct ProcessMemoryInfo
{
	size_t pageFaultCount;
	size_t workingSetSize;
	size_t quotaPagedPoolUsage;
	size_t quotaNonPagedPoolUsage;
	size_t pagefileUsage;
};

struct ProcessMemoryInfoEx : public ProcessMemoryInfo
{
	size_t peakWorkingSetSize;
	size_t quotaPeakPagedPoolUsage;
	size_t quotaPeakNonPagedPoolUsage;
	size_t peakPagefileUsage;
};

// https://learn.microsoft.com/en-us/windows/win32/psapi/collecting-memory-usage-information-for-a-process
ProcessMemoryInfo GetCurrentProcessMemoryInfo()
{
	ProcessMemoryInfo memInfo;

	PROCESS_MEMORY_COUNTERS pmc;
	if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
	{
		memInfo.pageFaultCount = pmc.PageFaultCount;
		memInfo.workingSetSize = pmc.WorkingSetSize;
		memInfo.quotaPagedPoolUsage = pmc.QuotaPagedPoolUsage;
		memInfo.quotaNonPagedPoolUsage = pmc.QuotaNonPagedPoolUsage;
		memInfo.pagefileUsage = pmc.PagefileUsage;
	}
	else
	{
		memset(&memInfo, 0, sizeof(memInfo));
	}

	return memInfo;
}

// https://learn.microsoft.com/en-us/windows/win32/psapi/collecting-memory-usage-information-for-a-process
ProcessMemoryInfoEx GetCurrentProcessMemoryInfoEx()
{
	ProcessMemoryInfoEx memInfo;

	PROCESS_MEMORY_COUNTERS pmc;
	if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
	{
		memInfo.pageFaultCount = pmc.PageFaultCount;
		memInfo.peakWorkingSetSize = pmc.PeakWorkingSetSize;
		memInfo.workingSetSize = pmc.WorkingSetSize;
		memInfo.quotaPeakPagedPoolUsage = pmc.QuotaPeakPagedPoolUsage;
		memInfo.quotaPagedPoolUsage = pmc.QuotaPagedPoolUsage;
		memInfo.quotaPeakNonPagedPoolUsage = pmc.QuotaPeakNonPagedPoolUsage;
		memInfo.quotaNonPagedPoolUsage = pmc.QuotaNonPagedPoolUsage;
		memInfo.pagefileUsage = pmc.PagefileUsage;
		memInfo.peakPagefileUsage = pmc.PeakPagefileUsage;
	}
	else
	{
		memset(&memInfo, 0, sizeof(memInfo));
	}

	return memInfo;
}

ProcessMemoryInfo operator+(ProcessMemoryInfo const & lhs, ProcessMemoryInfo const & rhs)
{
	ProcessMemoryInfo memInfo;
	memInfo.pageFaultCount         = lhs.pageFaultCount         + rhs.pageFaultCount;
	memInfo.workingSetSize         = lhs.workingSetSize         + rhs.workingSetSize;
	memInfo.quotaPagedPoolUsage    = lhs.quotaPagedPoolUsage    + rhs.quotaPagedPoolUsage;
	memInfo.quotaNonPagedPoolUsage = lhs.quotaNonPagedPoolUsage + rhs.quotaNonPagedPoolUsage;
	memInfo.pagefileUsage          = lhs.pagefileUsage          + rhs.pagefileUsage;
	return memInfo;
}

ProcessMemoryInfo operator-(ProcessMemoryInfo const & lhs, ProcessMemoryInfo const & rhs)
{
	ProcessMemoryInfo memInfo;
	memInfo.pageFaultCount         = lhs.pageFaultCount         - rhs.pageFaultCount;
	memInfo.workingSetSize         = lhs.workingSetSize         - rhs.workingSetSize;
	memInfo.quotaPagedPoolUsage    = lhs.quotaPagedPoolUsage    - rhs.quotaPagedPoolUsage;
	memInfo.quotaNonPagedPoolUsage = lhs.quotaNonPagedPoolUsage - rhs.quotaNonPagedPoolUsage;
	memInfo.pagefileUsage          = lhs.pagefileUsage          - rhs.pagefileUsage;
	return memInfo;
}

void DebugPrint(ProcessMemoryInfo const & memInfo)
{
	printf("\tPageFaultCount: %zd\n", memInfo.pageFaultCount);
	printf("\tWorkingSetSize: %s\n",  HumanReadableByteSizeStr(memInfo.workingSetSize).c_str());
	printf("\tQuotaPagedPoolUsage: %s\n", HumanReadableByteSizeStr(memInfo.quotaPagedPoolUsage).c_str());
	printf("\tQuotaNonPagedPoolUsage: %s\n", HumanReadableByteSizeStr(memInfo.quotaNonPagedPoolUsage).c_str());
	printf("\tPagefileUsage: %s\n", HumanReadableByteSizeStr(memInfo.pagefileUsage).c_str()); 
}

void DebugPrint(ProcessMemoryInfoEx const & memInfo)
{
	printf("\tPageFaultCount: %zd\n", memInfo.pageFaultCount);
	printf("\tPeakWorkingSetSize: %s\n", HumanReadableByteSizeStr(memInfo.peakWorkingSetSize).c_str());
	printf("\tWorkingSetSize: %s\n",  HumanReadableByteSizeStr(memInfo.workingSetSize).c_str());
	printf("\tQuotaPeakPagedPoolUsage: %s\n", HumanReadableByteSizeStr(memInfo.quotaPeakPagedPoolUsage).c_str());
	printf("\tQuotaPagedPoolUsage: %s\n", HumanReadableByteSizeStr(memInfo.quotaPagedPoolUsage).c_str());
	printf("\tQuotaPeakNonPagedPoolUsage: %s\n", HumanReadableByteSizeStr(memInfo.quotaPeakNonPagedPoolUsage).c_str());
	printf("\tQuotaNonPagedPoolUsage: %s\n", HumanReadableByteSizeStr(memInfo.quotaNonPagedPoolUsage).c_str());
	printf("\tPagefileUsage: %s\n", HumanReadableByteSizeStr(memInfo.pagefileUsage).c_str()); 
	printf("\tPeakPagefileUsage: %s\n", HumanReadableByteSizeStr(memInfo.peakPagefileUsage).c_str());
}


template <typename Func>
ProcessMemoryInfo MeasureFuncMemoryUsage(Func && f)
{
	const ProcessMemoryInfo prevProcessMemInfo = GetCurrentProcessMemoryInfo();
	f();
	return GetCurrentProcessMemoryInfo() - prevProcessMemInfo;
}


int BinaryFormatTests();
int MemoryMapTests();

int main()
{
	return BinaryFormatTests();
	// return MemoryMapTests();
}


#pragma region BinaryFormatTests

#include <cassert>
#include <vector>

#define K_USE_LZ4_COMPRESSION 1

#if K_USE_LZ4_COMPRESSION
#include <lz4.h>
#endif

// Inspired by https://github.com/raysan5/rres
struct BinDataHeader
{
	enum { MAX_CHUNKS = 1u << 16 };

	// File identifier
	union
	{
		unsigned char fourcc[4];
		uint32_t id;
	};

	uint16_t version; // File version: 100 for version 1.0
	uint16_t numChunks; // Number of data chunks in the file (in [1, 65536] -> store "actual number of chunks" - 1 in [0, 65535])
	uint32_t chunkUncompressedByteSize; // Byte size of an uncompressed chunk
	uint32_t reserved;
	uint64_t compressedByteSize;   // Compressed data byte size   (w/o BinDataHeader and BinDataChunkHeader(s))
	uint64_t uncompressedByteSize; // Uncompressed data byte size (w/o BinDataHeader and BinDataChunkHeader(s))
};

static_assert(sizeof(BinDataHeader) == 32, "Invalid BinDataHeader byte size");

struct BinDataChunkHeader
{
	uint64_t byteOffset;
	uint32_t compressedByteSize;
	uint32_t reserved;
};

static_assert(sizeof(BinDataChunkHeader) == 16, "Invalid BinDataChunkHeader byte size");


namespace details {

uint32_t GetChunkUncompressedByteSize(uint64_t srcByteSize)
{
	// Constraint: must generate at most BinDataHeader::MAX_CHUNKS chunks
	const uint64_t maxNumChunks = BinDataHeader::MAX_CHUNKS;
	uint32_t size = 4ull * 1024 * 1024; // 4Mo // TODO: tune this (make it dependenton srcByteSize?)
	uint64_t numChunks = (srcByteSize + size - 1) / size;
	while(numChunks > maxNumChunks)
	{
		size <<= 1;
		numChunks = (srcByteSize + size - 1) / size;
	}
	return size;
}

uint64_t GetCompressionByteSizeUpperBound(uint64_t srcByteSize)
{
	#if K_USE_LZ4_COMPRESSION
	return LZ4_compressBound((int)srcByteSize);
	#else
	return srcByteSize; // Dummy compression test
	#endif
}

}

uint64_t Compress(void * dst, void const * src, uint64_t srcByteSize, uint64_t dstMaxByteSize)
{
	#if K_USE_LZ4_COMPRESSION
	return LZ4_compress_default((char const *)src, (char *)dst, (int)srcByteSize, (int)dstMaxByteSize);
	#else
	// Dummy compression test function
	uint64_t const * iSrc = (uint64_t const *)src;
	uint32_t * iDst = (uint32_t *)dst;

	const uint64_t n = srcByteSize / sizeof(uint64_t);
	for(uint64_t i = 0; i < n; ++i)
	{
		iDst[i] = (uint32_t)iSrc[i];
	}

	return srcByteSize / 2;
	#endif
}

uint64_t Uncompress(void * dst, void const * src, uint64_t srcByteSize, uint64_t dstMaxByteSize)
{
	#if K_USE_LZ4_COMPRESSION
	return LZ4_decompress_safe((char const *)src, (char *)dst, srcByteSize, dstMaxByteSize);
	#else
	// Dummy decompression test function
	uint32_t const * iSrc = (uint32_t const *)src;
	uint64_t * iDst = (uint64_t *)dst;

	const uint64_t n = srcByteSize / sizeof(uint32_t);
	for(uint64_t i = 0; i < n; ++i)
	{
		iDst[i] = (uint64_t)iSrc[i];
	}

	return srcByteSize * 2;
	#endif
}

uint64_t GetCompressedBinaryDataByteSizeUpperBound(uint64_t uncompressedByteSize)
{
	const uint32_t chunkUncompressedByteSize = details::GetChunkUncompressedByteSize(uncompressedByteSize);
	const uint64_t numChunks = (uncompressedByteSize + chunkUncompressedByteSize - 1) / chunkUncompressedByteSize;

	const uint32_t chunkCompressedByteSizeUpperBound = details::GetCompressionByteSizeUpperBound(chunkUncompressedByteSize);

	return sizeof(BinDataHeader) + numChunks * (sizeof(BinDataChunkHeader) + chunkCompressedByteSizeUpperBound);
}

// \brief Compress binary data
// \param compressedData: destination buffer into which to write compressed data
// \param uncompressedData: source data to compress
// \param uncompressedDataByteSize: data to compress byte size
// \param compressedDataMaxByteSize: compressed data buffer byte size (should be at least equal to the value returned by GetCompressedBinaryDataByteSizeUpperBound)
// \return The size of the compressed binary data (<= value returned by GetCompressedBinaryDataByteSizeUpperBound) if successful, 0 otherwise.
uint64_t CompressedBinaryData(void * compressedData, void const * uncompressedData, uint64_t uncompressedDataByteSize, uint64_t compressedDataMaxByteSize)
{
	if(compressedDataMaxByteSize < GetCompressedBinaryDataByteSizeUpperBound(uncompressedDataByteSize))
		return 0;

	const uint32_t chunkUncompressedByteSize = details::GetChunkUncompressedByteSize(uncompressedDataByteSize);
	const uint64_t numChunks = (uncompressedDataByteSize + chunkUncompressedByteSize - 1) / chunkUncompressedByteSize;

	assert(numChunks > 0);
	assert(numChunks <= BinDataHeader::MAX_CHUNKS);

	const uint32_t chunkCompressedByteSizeUpperBound = details::GetCompressionByteSizeUpperBound(chunkUncompressedByteSize);

	const uint64_t compressedByteSizeUpperBound = sizeof(BinDataHeader) + numChunks * chunkCompressedByteSizeUpperBound;

	unsigned char * const baseCompressedData = (unsigned char *)compressedData;

	BinDataHeader * pHeader = (BinDataHeader*)baseCompressedData;
	pHeader->id = 0xC0DECAD0;
	pHeader->version = 100;
	pHeader->numChunks = static_cast<uint16_t>(numChunks-1);
	pHeader->chunkUncompressedByteSize = chunkUncompressedByteSize;
	pHeader->reserved = 0xDEADDEAD;
	pHeader->uncompressedByteSize = uncompressedDataByteSize;

	BinDataChunkHeader * pChunksHeaders = (BinDataChunkHeader *)(baseCompressedData + sizeof(BinDataHeader));
	unsigned char * pBaseDstChunksData = baseCompressedData + sizeof(BinDataHeader) + numChunks * sizeof(BinDataChunkHeader);

	unsigned char const * pBaseSrcData = (unsigned char const *)uncompressedData;
	unsigned char * pDstChunkData = pBaseDstChunksData;
	uint64_t chunkByteOffset = 0;
	for(uint64_t i = 0; i < numChunks; ++i)
	{
		pChunksHeaders[i].byteOffset = chunkByteOffset;
		unsigned char const * pSrcChunkData = pBaseSrcData + i * chunkUncompressedByteSize;
		const uint32_t srcChunkByteSize = i < numChunks-1 ? chunkUncompressedByteSize : uncompressedDataByteSize - i * chunkUncompressedByteSize;
		const uint64_t chunkCompressedByteSize = Compress(pDstChunkData, pSrcChunkData, srcChunkByteSize, chunkCompressedByteSizeUpperBound);
		assert(chunkCompressedByteSize < uint32_t(-1));
		pChunksHeaders[i].compressedByteSize = (uint32_t)chunkCompressedByteSize;
		pDstChunkData   += chunkCompressedByteSize;
		chunkByteOffset += chunkCompressedByteSize;
	}
	
	pHeader->compressedByteSize = chunkByteOffset;

	assert(pHeader->compressedByteSize == pChunksHeaders[numChunks-1].byteOffset + pChunksHeaders[numChunks-1].compressedByteSize);

	return sizeof(BinDataHeader) + pHeader->numChunks * sizeof(BinDataChunkHeader) + pHeader->compressedByteSize;

	// Note: could be compressed in parallel but would need compression at the end since each chunk would have allocated its upperbound
	#if 0
	unsigned char const * pBaseSrcData = (unsigned char const *)uncompressedData;
	parallel_for(uint64_t i = 0; i < numChunks; ++i)
	{
		unsigned char const * pSrcChunkData = pBaseSrcData + i * chunkUncompressedByteSize;
		unsigned char * pDstChunkData = pBaseDstChunksData + i * chunkCompressedByteSizeUpperBound;
		const uint32_t chunckUncompressedByteSize = i < numChunks-1 ? chunkUncompressedByteSize : uncompressedDataByteSize - i * chunkUncompressedByteSize;
		const uint64_t chunkCompressedByteSize = Compress(pSrcChunkData, pDstChunkData, chunckUncompressedByteSize, chunkCompressedByteSizeUpperBound);
		assert(chunkCompressedByteSize < uint32_t(-1));
		pChunksHeaders[i].compressedByteSize = (uint32_t)chunkCompressedByteSize;
	}

	// Setup chunks byte offsets (prefix scan)
	uint64_t chunkByteOffset = 0;
	for(uint64_t i = 0; i < numChunks; ++i)
	{
		pChunksHeaders[i].byteOffset = chunkByteOffset;
		chunkByteOffset += pChunksHeaders[i].compressedByteSize;
	}

	// Compact compressed chunks
	for(uint64_t i = 1; i < numChunks; ++i)
	{
		unsigned char const * pSrcChunkData = pBaseDstChunksData + i * chunkCompressedByteSizeUpperBound;
		unsigned char * pDstChunkData = pBaseDstChunksData + pChunksHeaders[i].byteOffset;
		memmove(pDstChunkData, pSrcChunkData, pChunksHeaders[i].compressedByteSize);
	}

	pHeader->compressedByteSize = chunkByteOffset;
	#endif
}

BinDataHeader const & GetBinDataHeader(void const * compressedData)
{
	return *(BinDataHeader const *)compressedData;
}

uint64_t GetUncompressedBinaryDataByteSize(void const * compressedData)
{
	return GetBinDataHeader(compressedData).uncompressedByteSize;
}

uint64_t GetCompressedBinaryDataNumChunks(void const * compressedData)
{
	return uint64_t(GetBinDataHeader(compressedData).numChunks) + 1;
}


using ChunkDecompressionFunc = void (*)(unsigned char const * pSrcChunkData, uint32_t srcChunkByteSize, uint64_t dstByteOffset, uint32_t dstChunkByteSize, void * userdata);

uint64_t UncompressedBinaryData(void const * compressedData, uint64_t compressedDataByteSize, ChunkDecompressionFunc chunkDecompressionFunc, void * userdata = nullptr)
{
	unsigned char const * pBaseCompressedData = (unsigned char const *)compressedData;

	BinDataHeader const * pHeader = (BinDataHeader const *)compressedData;

	const uint64_t dstUncompressedByteSize = pHeader->uncompressedByteSize;

	const uint64_t numChunks = uint64_t(pHeader->numChunks) + 1;
	const uint32_t chunkUncompressedByteSize = pHeader->chunkUncompressedByteSize;

	BinDataChunkHeader const * pChunksHeaders = (BinDataChunkHeader *)(pBaseCompressedData + sizeof(BinDataHeader));
	unsigned char const * pBaseCompressedChunksData = pBaseCompressedData + sizeof(BinDataHeader) + numChunks * sizeof(BinDataChunkHeader);

	for(uint64_t i = 0; i < numChunks; ++i)
	{
		unsigned char const * pSrcChunkData = pBaseCompressedChunksData + pChunksHeaders[i].byteOffset;
		const uint32_t srcChunkByteSize = pChunksHeaders[i].compressedByteSize;

		const uint64_t dstByteOffset = i * chunkUncompressedByteSize;
		const uint32_t dstChunkByteSize = i < numChunks-1 ? chunkUncompressedByteSize : dstUncompressedByteSize - i * chunkUncompressedByteSize;

		chunkDecompressionFunc(pSrcChunkData, srcChunkByteSize, dstByteOffset, dstChunkByteSize, userdata);
	}

	return dstUncompressedByteSize;
}

// \brief Uncompress binary data
// \param uncompressedData: destination buffer into which to write uncompressed data
// \param compressedData: source data to uncompress
// \param compressedDataByteSize: data to uncompress byte size
// \param uncompressedDataMaxByteSize: uncompressed data buffer byte size (should be at least equal to the value returned by GetUncompressedBinaryDataByteSize)
// \return The size of the uncompressed binary data (equal to value returned by GetUncompressedBinaryDataByteSize) if successful, 0 otherwise.
uint64_t UncompressedBinaryData(void * uncompressedData, void const * compressedData, uint64_t compressedDataByteSize, uint64_t uncompressedDataMaxByteSize)
{
	void * userdata = uncompressedData;
	const auto chunkDecompressionFunc = [](unsigned char const * pSrcChunkData, uint32_t srcChunkByteSize, uint64_t dstByteOffset, uint32_t dstChunkByteSize, void * userdata)
	{
		unsigned char * pDstChunkData = (unsigned char*)userdata + dstByteOffset;
		Uncompress(pDstChunkData, pSrcChunkData, srcChunkByteSize, dstChunkByteSize);
	};
	return UncompressedBinaryData(compressedData, compressedDataByteSize, chunkDecompressionFunc, userdata);
}

void PrintCompressedBinaryDataInfo(void * compressedData)
{
	BinDataHeader const & header = *(BinDataHeader const *)compressedData;
	printf("BinDataHeader:\n");
	printf("\tID: 0x%X\n", header.id);
	printf("\tVersion: %d\n", header.version);
	printf("\tUncompressed data byte size: %s\n", HumanReadableByteSizeStr(header.uncompressedByteSize).c_str());
	printf("\tCompressed data byte size: %s\n", HumanReadableByteSizeStr(header.compressedByteSize).c_str());
	printf("\tNum chunks: %llu\n", uint64_t(header.numChunks)+1);
	printf("\tChunk uncompressed byte size: %s\n", HumanReadableByteSizeStr(header.chunkUncompressedByteSize).c_str());
}

int BinaryFormatTests()
{
	#if ENABLE_PIX_RUNTIME
	HMODULE pixModule = PIXLoadLatestWinPixTimingCapturerLibrary();
	if(!pixModule)
	{
		DWORD errorCode = GetLastError();
		fprintf(stderr, "PIXLoadLatestWinPixTimingCapturerLibrary failed with error code: %d\n", errorCode);
		return 1;
	}

	wchar_t captureFilename[256];
	std::time_t time = std::time(nullptr);
	std::tm * calendarTime = std::localtime(&time);
	swprintf(captureFilename, sizeof(captureFilename)/sizeof(captureFilename[0]),
		L"BinaryFormatTests(%d-%02d-%02d.%02d-%02d-%02d).wpix",
		calendarTime->tm_year + 1900,
		calendarTime->tm_mon + 1,
		calendarTime->tm_mday,
		calendarTime->tm_hour,
		calendarTime->tm_min,
		calendarTime->tm_sec
	);

	PIXCaptureParameters params;
	params.TimingCaptureParameters.FileName = captureFilename;
	params.TimingCaptureParameters.CaptureCallstacks = true;
	params.TimingCaptureParameters.CaptureCpuSamples = 8000u; // must be 1000u, 4000, or 8000u. It's otherwise ignored.
	params.TimingCaptureParameters.CaptureFileIO = true;
	params.TimingCaptureParameters.CaptureVirtualAllocEvents = true;
	params.TimingCaptureParameters.CaptureHeapAllocEvents = true;
	params.TimingCaptureParameters.CaptureStorage = PIXCaptureParameters::PIXCaptureStorage::Memory;

	PIXBeginCapture(PIX_CAPTURE_TIMING, &params);
	#endif

	K_PIX_BEGIN_EVENT(0xFFFFFFFF, "BinaryFormatTests");

	K_PIX_BEGIN_EVENT(0xFFFF0000, "Init data");

	const uint64_t N = (1ull * 1024 * 1024 * 1024) / sizeof(uint64_t);
	const uint64_t srcDataByteSize = N * sizeof(uint64_t);
	uint64_t * srcData = (uint64_t *)malloc(srcDataByteSize);
	for(uint64_t i = 0; i < N; ++i)
		srcData[i] = i;

	K_PIX_END_EVENT();

	K_PIX_BEGIN_EVENT(0xFFFF8800, "Allocate compression buffer");

	const uint64_t compressedBinaryDataByteSizeUpperBound = GetCompressedBinaryDataByteSizeUpperBound(srcDataByteSize);
	void * compressionBuffer = malloc(compressedBinaryDataByteSizeUpperBound);
	
	K_PIX_END_EVENT();

	K_PIX_BEGIN_EVENT(0xFFFFFF00, "Compress data");

	auto start = std::chrono::high_resolution_clock::now();
	const uint64_t compressedBinaryDataByteSize = CompressedBinaryData(compressionBuffer, srcData, srcDataByteSize, compressedBinaryDataByteSizeUpperBound);
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	double duration_ms = duration.count() / 1000.0;
	printf("Compression time: %fms - [standard]\n", duration_ms);

	printf("Src byte size: %s\n", HumanReadableByteSizeStr(srcDataByteSize).c_str());
	printf("Dst byte size: %s\n", HumanReadableByteSizeStr(compressedBinaryDataByteSize).c_str());
	printf("Compression ratio: %f%%\n", double(compressedBinaryDataByteSize)/srcDataByteSize * 100.0);
	PrintCompressedBinaryDataInfo(compressionBuffer);
	
	K_PIX_END_EVENT();

	const uint64_t uncompressedBinaryDataByteSize = GetUncompressedBinaryDataByteSize(compressionBuffer);
	assert(uncompressedBinaryDataByteSize == srcDataByteSize);

	void * decompressionBuffer = malloc(uncompressedBinaryDataByteSize);

	// Standard
	K_PIX_BEGIN_EVENT(0xFF00FF00, "Uncompress data - standard");
	{
		memset(decompressionBuffer, 0xFF, uncompressedBinaryDataByteSize);

		start = std::chrono::high_resolution_clock::now();
		const uint64_t uncompressedBinaryDataByteSize2 = UncompressedBinaryData(decompressionBuffer, compressionBuffer, compressedBinaryDataByteSize, uncompressedBinaryDataByteSize);
		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		duration_ms = duration.count() / 1000.0;

		printf("Decompression time: %fms - [standard]\n", duration_ms);

		assert(uncompressedBinaryDataByteSize2 == uncompressedBinaryDataByteSize);
		assert(uncompressedBinaryDataByteSize2 == srcDataByteSize);

		uint64_t const * uncompressedData = (uint64_t const *)decompressionBuffer;
		int cmp = memcmp(srcData, uncompressedData, srcDataByteSize);
		if(cmp != 0)
		{
			fprintf(stderr, "Error while decompressing: uncompressedData != srcData [standard]\n");
		}
	}
	K_PIX_END_EVENT();

	// Parallel
	K_PIX_BEGIN_EVENT(0xFF0000FF, "Uncompress data - parallel");
	{
		memset(decompressionBuffer, 0xFF, uncompressedBinaryDataByteSize);

		struct DataChunkToProcess
		{
			unsigned char * dst;
			unsigned char const * src;
			uint32_t srcByteSize;
			uint32_t dstMaxByteSize;
		};

		// TODO: use SPMC queue?
		std::queue<DataChunkToProcess> dataChunksToProcess;

		std::mutex queueMutex;
		std::condition_variable queueCV;
		bool bReady = false;
		bool bKillThreads = false;

		const int numThreads = std::thread::hardware_concurrency() - 2;
		printf("Num threads: %d\n", numThreads);
		std::vector<std::thread> workers(numThreads);
		for(int i = 0; i < numThreads; ++i)
		{
			workers[i] = std::thread([&, i]()
			{
				while(true)
				{
					std::unique_lock<std::mutex> lock(queueMutex);
					queueCV.wait(lock, [&]{ return (bReady && !dataChunksToProcess.empty()) || bKillThreads; });

					if(bKillThreads)
						break;

					DataChunkToProcess chunk = dataChunksToProcess.front();
					dataChunksToProcess.pop();
					lock.unlock();

					K_PIX_BEGIN_EVENT(0xFF00FFFF, "Uncompress chunk - parallel");
					Uncompress(chunk.dst, chunk.src, chunk.srcByteSize, chunk.dstMaxByteSize);
					K_PIX_END_EVENT();
				}
			});
		}

		start = std::chrono::high_resolution_clock::now();

		const uint64_t numChunks = GetCompressedBinaryDataNumChunks(compressionBuffer);
		//dataChunksToProcess.reserve(numChunks);

		struct Userdata
		{
			void * pBaseDstChunkData;
			std::queue<DataChunkToProcess> & dataChunksToProcess;
		};

		Userdata userdata = { decompressionBuffer, dataChunksToProcess };
		const auto chunkDecompressionFunc = [](unsigned char const * pSrcChunkData, uint32_t srcChunkByteSize, uint64_t dstByteOffset, uint32_t dstChunkByteSize, void * pUserdata)
		{
			Userdata const & userdata = *(Userdata const *)pUserdata;
			unsigned char * pDstChunkData = (unsigned char*)userdata.pBaseDstChunkData + dstByteOffset;
			userdata.dataChunksToProcess.push({pDstChunkData, pSrcChunkData, srcChunkByteSize, dstChunkByteSize});
		};

		K_PIX_BEGIN_EVENT(0xFF00FFFF, "Prepare chunks - parallel");
		const uint64_t uncompressedBinaryDataByteSize2 = UncompressedBinaryData(compressionBuffer, compressedBinaryDataByteSize, chunkDecompressionFunc, &userdata);
		K_PIX_END_EVENT();

		#if 0
		while(!dataChunksToProcess.empty())
		{
			DataChunkToProcess const & chunk = dataChunksToProcess.front();
			Uncompress(chunk.dst, chunk.src, chunk.srcByteSize, chunk.dstMaxByteSize);
			dataChunksToProcess.pop();
		}
		#else
		// Notify workers to start working
		{
			std::unique_lock<std::mutex> lock(queueMutex);
			bReady = true;
		}
		queueCV.notify_all();
		#endif

		// Waits for workers to finish
		while(true)
		{
			bool bProcessingDone;
			{
				std::unique_lock<std::mutex> lock(queueMutex);
				if(bProcessingDone = dataChunksToProcess.empty())
				{
					bKillThreads = true;
					queueCV.notify_all();
				}
			}
			
			if(bProcessingDone)
				break;

			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}

		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		duration_ms = duration.count() / 1000.0;
		printf("Decompression time: %fms - [parallel]\n", duration_ms);

		for(std::thread & w : workers)
			w.join();

		assert(uncompressedBinaryDataByteSize2 == uncompressedBinaryDataByteSize);
		assert(uncompressedBinaryDataByteSize2 == srcDataByteSize);

		// TODO: try to notify threads while main thread queue jobs

		uint64_t const * uncompressedData = (uint64_t const *)decompressionBuffer;
		int cmp = memcmp(srcData, uncompressedData, srcDataByteSize);
		if(cmp != 0)
		{
			fprintf(stderr, "Error while decompressing: uncompressedData != srcData [parallel]\n");
		}
	}
	K_PIX_END_EVENT();

	// Memcpy
	{
		const auto start = std::chrono::high_resolution_clock::now();
		memcpy(decompressionBuffer, srcData, srcDataByteSize);
		const auto stop = std::chrono::high_resolution_clock::now();
		const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		const double duration_ms = duration.count() / 1000.0;
		printf("Decompression time: %fms - [memcpy]\n", duration_ms);
	}

	free(srcData);
	free(compressionBuffer);
	free(decompressionBuffer);

	K_PIX_END_EVENT();

	#if ENABLE_PIX_RUNTIME
	PIXEndCapture(false);

	if(!FreeLibrary(pixModule))
		fprintf(stderr, "Failed to free PIX library\n");
	#endif

	return 0;
}

#pragma endregion

#pragma region MemoryMapTests

// mmap:
//  pros: 
//    + easier multithreading -> map file, access buffer from multiple threads
//    + caching/swapping automatically handle by OS
//  cons:
//    - potentially higher memory usage
//    - higher TLB miss because mmap does not wire all the page (unless w/ MAP_POPULATE)

// https://yarchive.net/comp/mmap.html#6

int MemoryMapTests()
{
	char const * filename = "data.bin";

	// Write a 1GB binary file
	if(0)
	{
		const uint64_t baseFileByteSize = 1024 * 1024 * 1024;

		// https://learn.microsoft.com/en-us/windows/win32/memory/creating-a-view-within-a-file
		// https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-createfilea
		HANDLE hFile = CreateFile(filename,
			GENERIC_READ | GENERIC_WRITE,
			0,
			NULL,
			CREATE_ALWAYS,
			FILE_ATTRIBUTE_NORMAL,
			NULL
		);

		if(hFile == INVALID_HANDLE_VALUE)
		{
			fprintf(stderr, "Failed to create file: '%s'", filename);
			return 1;
		}

		// Get the system allocation granularity
		SYSTEM_INFO sysInfo;
		GetSystemInfo(&sysInfo);
		DWORD dwSysGran = sysInfo.dwAllocationGranularity;

		const uint64_t baseFileMapStart = 0;

		// To calculate where to start the file mapping, round down the
		// offset of the data into the file to the nearest multiple of the
		// system allocation granularity.
		const uint64_t fileMapStart = (baseFileMapStart / dwSysGran) * dwSysGran;
		DWORD fileMapStartHi = (DWORD)(fileMapStart >> 32);
		DWORD fileMapStartLo = (DWORD)(fileMapStart & 0xFFFFFFFF);
		fprintf(stderr, "The file map view starts at %llu bytes into the file.\n", fileMapStart);

		// Calculate the size of the file mapping view.
		const uint64_t mapViewSize = (baseFileMapStart % dwSysGran) + baseFileByteSize;
		const DWORD mapViewSizeHi = (DWORD)(mapViewSize >> 32);
		const DWORD mapViewSizeLo = (DWORD)(mapViewSize & 0xFFFFFFFF);
		fprintf(stderr, "The file map view is %llu bytes large.\n", mapViewSize);

		// How large will the file mapping object be?
		const size_t fileMapSize = baseFileMapStart + baseFileByteSize;
		const DWORD fileMapSizeHi = (DWORD)(fileMapSize >> 32);
		const DWORD fileMapSizeLo = (DWORD)(fileMapSize & 0xFFFFFFFF);
		fprintf(stderr, "The file mapping object is %llu bytes large.\n", fileMapSize);

		// The data of interest isn't at the beginning of the
		// view, so determine how far into the view to set the pointer.
		const uint64_t iViewDelta = baseFileMapStart - fileMapStart;
		fprintf(stderr, "The data is %llu bytes into the view.\n", iViewDelta);
		
		#if 0
		// Now write a file with data suitable for experimentation. This
		// provides unique int (4-byte) offsets in the file for easy visual
		// inspection. Note that this code does not check for storage
		// medium overflow or other errors, which production code should
		// do. Because an int is 4 bytes, the value at the pointer to the
		// data should be one quarter of the desired offset into the file
		DWORD dBytesWritten = 0;
		for(int i = 0; i < (int)dwSysGran; ++i)
		{
			WriteFile(hFile, &i, sizeof(i), &dBytesWritten, NULL);
		}

		// Verify that the correct file size was written.
		DWORD fileSizeHi = 0;
		DWORD fileSizeLo = GetFileSize(hFile, &fileSizeHi);
		uint64_t fileSize = (uint64_t(fileSizeHi) << 32) | fileSizeLo;
		fprintf(stderr, "hFile size: %llu\n", fileSize);
		#endif

		// Create a file mapping object for the file
		// Note that it is a good idea to ensure the file size is not zero
		HANDLE hMapFile = CreateFileMapping(
			hFile,          // current file handle
			NULL,           // default security
			PAGE_READWRITE, // read/write permission
			mapViewSizeHi,  // size of mapping object, high
			mapViewSizeLo,  // size of mapping object, low
			NULL            // name of mapping object
		);

		if(hMapFile == NULL)
		{
			fprintf(stderr, "Failed to create mapping of file: '%s'", filename);
			return 2;
		}

		// Map the view and test the results.
		LPVOID lpMapAddress = MapViewOfFile(
			hMapFile,            // handle to mapping object
			FILE_MAP_ALL_ACCESS, // read/write
			fileMapStartHi,      // high-order 32 bits of file offset
			fileMapStartLo,      // low-order 32 bits of file offset
			mapViewSize          // number of bytes to map
		);

		if(lpMapAddress == NULL)
		{
			fprintf(stderr, "lpMapAddress is NULL: last error: %d\n", GetLastError());
			return 3;
		}

		// Verify that the correct file size was written.
		DWORD fileSizeHi = 0;
		DWORD fileSizeLo = GetFileSize(hFile, &fileSizeHi);
		uint64_t fileSize = (uint64_t(fileSizeHi) << 32) | fileSizeLo;
		fprintf(stderr, "hFile size: %lluB\n", fileSize);

		// Calculate the pointer to the data.
		char * pData = (char *)lpMapAddress + iViewDelta;

		// Extract the data, an int. Cast the pointer pData from a "pointer
		// to char" to a "pointer to int" to get the whole thing
		int * iData = (int *)pData;

		const uint64_t numInts = baseFileByteSize / sizeof(float);
		for(int i = 0; i < (int)numInts; ++i)
		{
			iData[i] = i;
		}

		if(!UnmapViewOfFile(lpMapAddress))
		{
			fprintf(stderr, "UnmapViewOfFile error: %d\n", GetLastError());
		}

		if(!CloseHandle(hMapFile))
		{
			fprintf(stderr, "CloseHandle(hMapFile) error: %d\n", GetLastError());
		}

		if(!CloseHandle(hFile))
		{
			fprintf(stderr, "CloseHandle(hFile) error: %d\n", GetLastError());
		}

		return 0;
	}

	// Read file
	if(1)
	{
		#if ENABLE_PIX_RUNTIME
		HMODULE pixModule = PIXLoadLatestWinPixTimingCapturerLibrary();
		if(!pixModule)
		{
			DWORD errorCode = GetLastError();
			fprintf(stderr, "PIXLoadLatestWinPixTimingCapturerLibrary failed with error code: %d\n", errorCode);
			return 1;
		}

		wchar_t captureFilename[256];
		std::time_t time = std::time(nullptr);
		std::tm * calendarTime = std::localtime(&time);
		swprintf(captureFilename, sizeof(captureFilename)/sizeof(captureFilename[0]),
			L"MmapTest(%d-%02d-%02d.%02d-%02d-%02d).wpix",
			calendarTime->tm_year + 1900,
			calendarTime->tm_mon + 1,
			calendarTime->tm_mday,
			calendarTime->tm_hour,
			calendarTime->tm_min,
			calendarTime->tm_sec
		);

		PIXCaptureParameters params;
		params.TimingCaptureParameters.FileName = captureFilename;
		params.TimingCaptureParameters.CaptureCallstacks = true;
		params.TimingCaptureParameters.CaptureCpuSamples = 8000u; // must be 1000u, 4000, or 8000u. It's otherwise ignored.
		params.TimingCaptureParameters.CaptureFileIO = true;
		params.TimingCaptureParameters.CaptureVirtualAllocEvents = true;
		params.TimingCaptureParameters.CaptureHeapAllocEvents = true;
		params.TimingCaptureParameters.CaptureStorage = PIXCaptureParameters::PIXCaptureStorage::Memory;

		PIXBeginCapture(PIX_CAPTURE_TIMING, &params);
		#endif

		PROCESS_MEMORY_COUNTERS pmc;

		// https://learn.microsoft.com/en-us/windows/win32/psapi/collecting-memory-usage-information-for-a-process
		if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
		{
			printf("\tPageFaultCount: %d\n", pmc.PageFaultCount);
			printf("\tPeakWorkingSetSize: %zd\n", pmc.PeakWorkingSetSize);
			printf("\tPeakWorkingSetSize: %s\n", HumanReadableByteSizeStr(pmc.PeakWorkingSetSize).c_str());
			printf("\tWorkingSetSize: %zd\n", pmc.WorkingSetSize);
			printf("\tWorkingSetSize: %s\n",  HumanReadableByteSizeStr(pmc.WorkingSetSize).c_str());
			printf("\tQuotaPeakPagedPoolUsage: %zd\n", pmc.QuotaPeakPagedPoolUsage);
			printf("\tQuotaPagedPoolUsage: %zd\n", pmc.QuotaPagedPoolUsage);
			printf("\tQuotaPeakNonPagedPoolUsage: %zd\n", pmc.QuotaPeakNonPagedPoolUsage);
			printf("\tQuotaNonPagedPoolUsage: %zd\n", pmc.QuotaNonPagedPoolUsage);
			printf("\tPagefileUsage: %zd\n", pmc.PagefileUsage); 
			printf("\tPagefileUsage: %s\n", HumanReadableByteSizeStr(pmc.PagefileUsage).c_str()); 
			printf("\tPeakPagefileUsage: %zd\n", pmc.PeakPagefileUsage);
			printf("\tPeakPagefileUsage: %s\n", HumanReadableByteSizeStr(pmc.PeakPagefileUsage).c_str());
		}

		#if 0
		K_PIX_BEGIN_EVENT(0xFFFFFFFF, "main - fread");
		{
			FILE * file = fopen(filename, "rb");

			fseek(file, 0L, SEEK_END);
			size_t size = ftell(file);
			fseek(file, 0L, SEEK_SET);

			void * data = malloc(size);
		
			printf("\n");
		
			if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
				printf("\tPageFaultCount: %d\n", pmc.PageFaultCount);
			K_PIX_BEGIN_EVENT(0xFFFF00FF, "fread");
			fread(data, 1, size, file);
			K_PIX_END_EVENT();

			if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
				printf("\tPageFaultCount: %d\n", pmc.PageFaultCount);

			int * iData = (int*)data;
			size_t sum = 0;
			for(size_t i = 0; i < size/sizeof(int); ++i)
			{
				sum += iData[i];
			}

			if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
				printf("\tPageFaultCount: %d\n", pmc.PageFaultCount);

			printf("sum: %zd\n", sum);

			free(data);
		
			fclose(file);
		}
		K_PIX_END_EVENT();
		#endif

		#if 1
		K_PIX_BEGIN_EVENT(0xFFFFFFFF, "main - mmap");
		{
			// https://learn.microsoft.com/en-us/windows/win32/memory/creating-a-view-within-a-file
			// https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-createfilea

			DWORD dwFlagsAndAttributes = FILE_ATTRIBUTE_NORMAL;
			dwFlagsAndAttributes = FILE_FLAG_SEQUENTIAL_SCAN;
			//dwFlagsAndAttributes = FILE_FLAG_RANDOM_ACCESS;

			HANDLE hFile = CreateFile(filename,
				GENERIC_READ,
				0,
				NULL,
				OPEN_EXISTING,
				dwFlagsAndAttributes,
				NULL
			);

			if(hFile == INVALID_HANDLE_VALUE)
			{
				fprintf(stderr, "Failed to create file: '%s'", filename);
				return 1;
			}

			// Create a file mapping object for the file
			// Note that it is a good idea to ensure the file size is not zero
			HANDLE hMapFile = CreateFileMapping(
				hFile,          // current file handle
				NULL,           // default security
				PAGE_READONLY,  // read/write permission
				0,  // size of mapping object, high
				0,  // size of mapping object, low
				NULL            // name of mapping object
			);

			if(hMapFile == NULL)
			{
				fprintf(stderr, "Failed to create mapping of file: '%s'", filename);
				return 2;
			}


			printf("\n");
		
			if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
			{
				printf("\tPageFaultCount: %d\n", pmc.PageFaultCount);
				printf("\tPeakWorkingSetSize: %s\n", HumanReadableByteSizeStr(pmc.PeakWorkingSetSize).c_str());
				printf("\tWorkingSetSize: %s\n",  HumanReadableByteSizeStr(pmc.WorkingSetSize).c_str());
			}

			// Map the view and test the results.
			// https://learn.microsoft.com/en-us/windows/win32/api/memoryapi/nf-memoryapi-mapviewoffile
			LPVOID lpMapAddress = MapViewOfFile(
				hMapFile,            // handle to mapping object
				FILE_MAP_READ, // read/write
				0,      // high-order 32 bits of file offset
				0,      // low-order 32 bits of file offset
				0          // number of bytes to map
			);

			if(lpMapAddress == NULL)
			{
				fprintf(stderr, "lpMapAddress is NULL: last error: %d\n", GetLastError());
				ErrorExit("MapViewOfFile");
				return 3;
			}

			// Verify that the correct file size was written.
			DWORD fileSizeHi = 0;
			DWORD fileSizeLo = GetFileSize(hFile, &fileSizeHi);
			uint64_t fileSize = (uint64_t(fileSizeHi) << 32) | fileSizeLo;
			fprintf(stderr, "hFile size: %lluB\n", fileSize);

			char * pData = (char *)lpMapAddress;

			// Extract the data, an int. Cast the pointer pData from a "pointer
			// to char" to a "pointer to int" to get the whole thing
			int * iData = (int *)pData;

			if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
			{
				printf("\tPageFaultCount: %d\n", pmc.PageFaultCount);
				printf("\tPeakWorkingSetSize: %s\n", HumanReadableByteSizeStr(pmc.PeakWorkingSetSize).c_str());
				printf("\tWorkingSetSize: %s\n",  HumanReadableByteSizeStr(pmc.WorkingSetSize).c_str());
			}

			size_t prevPageFaultCount = 0;
			if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
			{
				printf("\tPageFaultCount: %d\n", pmc.PageFaultCount);
				prevPageFaultCount = pmc.PageFaultCount;
			}


			// https://stackoverflow.com/questions/38381233/find-next-aligned-memory-address
			const auto Align = [](size_t alignment, void * ptr) -> void *
			{
				const std::uintptr_t intptr = reinterpret_cast<std::uintptr_t>(ptr);
				const std::uintptr_t aligned = (intptr - 1u + alignment) & -alignment;
				return reinterpret_cast<void *>(aligned);
			};

			const auto FindMaxAlignment = [&Align](void * ptr) -> size_t
			{
				const size_t kMaxAlignment = 1u << 16;
				size_t alignment = kMaxAlignment;
				do
				{
					void * alignedPtr = Align(alignment, ptr);
					if(alignedPtr == ptr)
						break;
				} while(alignment >>= 1);
				return alignment;
			};

			const size_t maxAlignment = FindMaxAlignment(iData);
			printf("0x%p | max alignment: %s\n", iData, HumanReadableByteSizeStr(maxAlignment).c_str());

			size_t sum = 0;
			ProcessMemoryInfo sumMemoryUsage = MeasureFuncMemoryUsage([&]()
			{
				for(size_t i = 0; i < fileSize/sizeof(int); ++i)
				{
					sum += iData[i];
				}
			});

			printf("\nmmap - sum\n");
			DebugPrint(sumMemoryUsage);

			// https://learn.microsoft.com/en-us/windows/win32/api/sysinfoapi/ns-sysinfoapi-system_info
			SYSTEM_INFO sysInfo;
			GetSystemInfo(&sysInfo);

			// The page size and the granularity of page protection and commitment.
			// This is the page size used by the VirtualAlloc function.
			const size_t pageSize = sysInfo.dwPageSize;

			// The granularity for the starting address at which virtual memory can be allocated.
			const size_t allocationGranularity = sysInfo.dwAllocationGranularity;

			size_t pageFaultSize = sumMemoryUsage.pageFaultCount * pageSize;
			printf("pageSize: %s\n", HumanReadableByteSizeStr(pageSize).c_str());
			printf("allocationGranularity: %s\n", HumanReadableByteSizeStr(allocationGranularity).c_str());
			printf("pageFaultSize: %s\n", HumanReadableByteSizeStr(pageFaultSize).c_str());
			printf("theorical page fault count: %zd\n", sumMemoryUsage.workingSetSize/pageSize);

			printf("\nsum: %zd\n", sum);

			if(!UnmapViewOfFile(lpMapAddress))
			{
				fprintf(stderr, "UnmapViewOfFile error: %d\n", GetLastError());
			}

			if(!CloseHandle(hMapFile))
			{
				fprintf(stderr, "CloseHandle(hMapFile) error: %d\n", GetLastError());
			}

			if(!CloseHandle(hFile))
			{
				fprintf(stderr, "CloseHandle(hFile) error: %d\n", GetLastError());
			}
		}
		K_PIX_END_EVENT();
		#endif

		// https://stackoverflow.com/questions/17925051/fast-textfile-reading-in-c
		// "Reading in 16kiB chunks reuses the same 4 pages of address space in your process.
		// You won't have TLB misses, and 16kiB is smaller than L1 cache.
		// The memcpy from page-cache (inside read(2)) goes very fast, and the memchr only touches memory that's hot in L1.
		// The mmap version has to fault each page, because mmap doesn't wire all the pages
		// (unless you use MAP_POPULATE, but that won't work well when file size is a large fraction of RAM size)"
		{
			#if 1
			FILE * file = fopen(filename, "rb");

			fseek(file, 0L, SEEK_END);
			size_t size = ftell(file);
			fseek(file, 0L, SEEK_SET);
			#else
			// https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-createfilea
			DWORD dwFlagsAndAttributes = FILE_ATTRIBUTE_NORMAL;
			dwFlagsAndAttributes = FILE_FLAG_SEQUENTIAL_SCAN;
			HANDLE hFile = CreateFile(filename,
				GENERIC_READ,
				0,
				NULL,
				OPEN_EXISTING,
				dwFlagsAndAttributes,
				NULL
			);
			#endif

			constexpr size_t kBufferSize = 16 * 1024; // 16kB // TODO: can we query L1 cache size?

			size_t prevPageFaultCount = 0;
			if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
			{
				printf("\tPageFaultCount: %d\n", pmc.PageFaultCount);
				prevPageFaultCount = pmc.PageFaultCount;
			}

			size_t sum = 0;
			char buffer[kBufferSize + 1];
			while(size_t bytes_read = fread(buffer, 1, kBufferSize, file))
			{
				if(!bytes_read)
				{
					if (feof(file))
						printf("Error reading file: unexpected end of file\n");
					else if (ferror(file))
						perror("Error reading file");
					break;
				}

				int * iData = (int*)buffer;
				for(size_t i = 0; i < bytes_read/sizeof(int); ++i)
				{
					sum += iData[i];
				}
			}

			printf("Sum: %zd\n", sum);

			size_t currPageFaultCount = 0;
			if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
			{
				printf("\tPageFaultCount: %d\n", pmc.PageFaultCount);
				currPageFaultCount = pmc.PageFaultCount;
			}

			size_t pageFaultCount = currPageFaultCount - prevPageFaultCount;
			printf("pageFaultCount: %zd\n", pageFaultCount);

			fclose(file);
		}

		#if 0
		// Get the list of process identifiers.
		DWORD aProcesses[1024], cbNeeded, cProcesses;
		unsigned int i;

		if ( !EnumProcesses( aProcesses, sizeof(aProcesses), &cbNeeded ) )
		{
			return 1;
		}

		// Calculate how many process identifiers were returned.
		cProcesses = cbNeeded / sizeof(DWORD);

		// Print the memory usage for each process
		for ( i = 0; i < cProcesses; i++ )
		{
			PrintMemoryInfo( aProcesses[i] );
		}
		#endif

		printf("\n");

		if(GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
		{
			printf("\tPageFaultCount: %d\n", pmc.PageFaultCount);
			printf("\tPeakWorkingSetSize: %zd\n", pmc.PeakWorkingSetSize);
			printf("\tPeakWorkingSetSize: %s\n", HumanReadableByteSizeStr(pmc.PeakWorkingSetSize).c_str());
			printf("\tWorkingSetSize: %zd\n", pmc.WorkingSetSize);
			printf("\tWorkingSetSize: %s\n",  HumanReadableByteSizeStr(pmc.WorkingSetSize).c_str());
			printf("\tQuotaPeakPagedPoolUsage: %zd\n", pmc.QuotaPeakPagedPoolUsage);
			printf("\tQuotaPagedPoolUsage: %zd\n", pmc.QuotaPagedPoolUsage);
			printf("\tQuotaPeakNonPagedPoolUsage: %zd\n", pmc.QuotaPeakNonPagedPoolUsage);
			printf("\tQuotaNonPagedPoolUsage: %zd\n", pmc.QuotaNonPagedPoolUsage);
			printf("\tPagefileUsage: %zd\n", pmc.PagefileUsage); 
			printf("\tPagefileUsage: %s\n", HumanReadableByteSizeStr(pmc.PagefileUsage).c_str()); 
			printf("\tPeakPagefileUsage: %zd\n", pmc.PeakPagefileUsage);
			printf("\tPeakPagefileUsage: %s\n", HumanReadableByteSizeStr(pmc.PeakPagefileUsage).c_str());
		}

		#if ENABLE_PIX_RUNTIME
		PIXEndCapture(false);

		if(!FreeLibrary(pixModule))
			fprintf(stderr, "Failed to free PIX library\n");
		#endif
	}
	


	#if 0
	CUdevice cuDevice;
	CUcontext context;
	CUDA_SAFE_CALL(cuInit(0));
	CUDA_SAFE_CALL(cuDeviceGet(&cuDevice, 0));
	CUDA_SAFE_CALL(cuCtxCreate(&context, 0, cuDevice));

	std::mutex mutex;
	std::condition_variable cv;
	bool bReady = false;

	const int numThreads = std::thread::hardware_concurrency() - 2;
	std::vector<std::thread> workers(numThreads);
	for(int i = 0; i < numThreads; ++i)
	{
		workers[i] = std::thread([&, i]()
		{
			CUDA_SAFE_CALL(cuCtxSetCurrent(context));

			std::unique_lock<std::mutex> lock(mutex);
			cv.wait(lock, [&]{ return bReady; });
			lock.unlock();
		});
	}

	#if ENABLE_PIX_RUNTIME
	HMODULE pixModule = PIXLoadLatestWinPixTimingCapturerLibrary();
	if(!pixModule)
	{
		DWORD errorCode = GetLastError();
		fprintf(stderr, "PIXLoadLatestWinPixTimingCapturerLibrary failed with error code: %d\n", errorCode);
		return 1;
	}

	wchar_t captureFilename[256];
	std::time_t time = std::time(nullptr);
	std::tm * calendarTime = std::localtime(&time);
	swprintf(captureFilename, sizeof(captureFilename)/sizeof(captureFilename[0]),
		L"MmapTest(%d-%02d-%02d.%02d-%02d-%02d).wpix",
		calendarTime->tm_year + 1900,
		calendarTime->tm_mon + 1,
		calendarTime->tm_mday,
		calendarTime->tm_hour,
		calendarTime->tm_min,
		calendarTime->tm_sec
	);

	PIXCaptureParameters params;
	params.TimingCaptureParameters.FileName = captureFilename;
	params.TimingCaptureParameters.CaptureCallstacks = true;
	params.TimingCaptureParameters.CaptureCpuSamples = 8000u; // must be 1000u, 4000, or 8000u. It's otherwise ignored.
	params.TimingCaptureParameters.CaptureFileIO = true;
	params.TimingCaptureParameters.CaptureVirtualAllocEvents = true;
	params.TimingCaptureParameters.CaptureHeapAllocEvents = true;
	params.TimingCaptureParameters.CaptureStorage = PIXCaptureParameters::PIXCaptureStorage::Memory;

	PIXBeginCapture(PIX_CAPTURE_TIMING, &params);
	#endif

	K_PIX_BEGIN_EVENT(0xFFFFFFFF, "main");

	{
		std::lock_guard<std::mutex> lock(mutex);
		bReady = true;
		cv.notify_all();
	}

	for(std::thread & w : workers)
		w.join();

	K_PIX_END_EVENT();

	#if ENABLE_PIX_RUNTIME
	PIXEndCapture(false);

	if(!FreeLibrary(pixModule))
		fprintf(stderr, "Failed to free PIX library\n");
	#endif

	CUDA_SAFE_CALL(cuCtxDestroy(context));
	#endif

	return 0;
}

#pragma endregion