#include <cuda.h>
#include <cuda_runtime.h>
#include <nvtx3/nvToolsExt.h>

#include <stdio.h>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

#include <Windows.h>
#include <psapi.h>

#undef ENABLE_PIX_RUNTIME
#define ENABLE_PIX_RUNTIME 0

#undef ENABLE_NVTX
#define ENABLE_NVTX 1

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

#if ENABLE_NVTX
#define K_NVTX_RANGE_PUSH(colorARGB, str) \
do \
{ \
	nvtxEventAttributes_t eventAttrib = {0}; \
	eventAttrib.version = NVTX_VERSION; \
	eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE; \
	eventAttrib.colorType = NVTX_COLOR_ARGB; \
	eventAttrib.color = colorARGB; \
	eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII; \
	eventAttrib.message.ascii = str; \
	nvtxRangePushEx(&eventAttrib); \
} while(0)

#define K_NVTX_RANGE_POP() nvtxRangePop()
#else
#define K_NVTX_RANGE_PUSH(colorARGB, str) K_NOOP()
#define K_NVTX_RANGE_POP()                K_NOOP()
#endif

#define K_PUSH_PROFILING_MARKER(colorARGB, str) K_PIX_BEGIN_EVENT(colorARGB, str); K_NVTX_RANGE_PUSH(colorARGB, str)
#define K_POP_PROFILING_MARKER()                K_PIX_END_EVENT(); K_NVTX_RANGE_POP()

#define CUDA_DRIVER_CHECK_CALL(x)                           \
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

#define CUDA_CHECK_CALL(call)                       \
  do {                                              \
	cudaError_t err = call;                         \
	if(err != cudaSuccess) {                        \
	  fprintf(stderr, "CUDA error at %s %d: %s\n",  \
		 __FILE__, __LINE__,                        \
		 cudaGetErrorString(err));                  \
	  exit(EXIT_FAILURE);                           \
	}                                               \
  } while (0)


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

template <typename Func>
double MeasureFuncTime(Func && f)
{
	const auto start = std::chrono::high_resolution_clock::now();
	f();
	const auto stop = std::chrono::high_resolution_clock::now();
	const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	const double duration_ms = duration.count() / 1000.0;
	return duration_ms;
}

template <typename Func>
double MeasureAndPrintFuncTime(char const * label, Func && f)
{
	const double duration_ms = MeasureFuncTime(std::forward<Func>(f));
	printf("%s: %fms\n", label, duration_ms);
	return duration_ms;
}


// https://github.com/kirksaunders/barrier/blob/master/barrier.hpp
// https://stackoverflow.com/questions/24465533/implementing-boostbarrier-in-c11/
// https://codereview.stackexchange.com/questions/243519/barrier-implementation-in-c
// https://stackoverflow.com/questions/17101922/do-i-have-to-acquire-lock-before-calling-condition-variable-notify-one
// https://stackoverflow.com/questions/26190853/why-does-the-boost-library-use-an-m-generation-variable-in-its-implementation-of
class Barrier
{
	public:
		// Construct barrier for use with num threads.
		Barrier(std::size_t num)
			: num_threads(num),
				wait_count(0),
				instance(0),
				mut(),
				cv()
		{

		}

		// disable copying of barrier
		Barrier(const Barrier&) = delete;
		Barrier& operator =(const Barrier&) = delete;

		// This function blocks the calling thread until
		// all threads (specified by num_threads) have
		// called it. Blocking is achieved using a
		// call to condition_variable.wait().
		void wait()
		{
			std::unique_lock<std::mutex> lock(mut); // acquire lock
			std::size_t inst = instance; // store current instance for comparison
											// in predicate

			if(++wait_count == num_threads) // all threads reached barrier
			{
				wait_count = 0; // reset wait_count
				instance++; // increment instance for next use of barrier and to
							// pass condition variable predicate
				cv.notify_all();
			}
			else // not all threads have reached barrier
			{
				cv.wait(lock, [this, &inst]() { return instance != inst; });
				// NOTE: The predicate lambda here protects against spurious
				//       wakeups of the thread. As long as this->instance is
				//       equal to inst, the thread will not wake.
				//       this->instance will only increment when all threads
				//       have reached the barrier and are ready to be unblocked.
			}
		}
	private:
		std::size_t num_threads; // number of threads using barrier
		std::size_t wait_count; // counter to keep track of waiting threads
		std::size_t instance; // counter to keep track of barrier use count
		std::mutex mut; // mutex used to protect resources
		std::condition_variable cv; // condition variable used to block threads
};


int PinnedHostMemoryTests();
int BinaryFormatTests();
int MemoryMapTests();

int main()
{
	//return PinnedHostMemoryTests();
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
	enum
	{
		MAX_CHUNKS = 1u << 16,
		MIN_CHUNKS_UNCOMPRESSED_BYTE_SIZE = 4096, // 4Ko
		MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE = 16 * 1024 * 1024 // 64Mo -> up to MAX_CHUNKS * MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE bytes = 65536 * 16Mo = 1To
	};

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
	uint32_t size = 4ull * 1024 * 1024; // 4Mo // TODO: tune this (make it dependent on srcByteSize?)
	uint64_t numChunks = (srcByteSize + size - 1) / size;
	while(numChunks > maxNumChunks)
	{
		size <<= 1;
		numChunks = (srcByteSize + size - 1) / size;
	}
	assert(size >= BinDataHeader::MIN_CHUNKS_UNCOMPRESSED_BYTE_SIZE);
	assert(size <= BinDataHeader::MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE);
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
	#if 1
		uint64_t const * iSrc = (uint64_t const *)src;
		uint32_t * iDst = (uint32_t *)dst;

		const uint64_t n = srcByteSize / sizeof(uint64_t);
		for(uint64_t i = 0; i < n; ++i)
		{
			iDst[i] = (uint32_t)iSrc[i];
		}

		return srcByteSize / 2;
	#else
		uint64_t const * iSrc = (uint64_t const *)src;
		uint64_t * iDst = (uint64_t *)dst;

		const uint64_t n = srcByteSize / sizeof(uint64_t);
		for(uint64_t i = 0; i < n; ++i)
		{
			iDst[i] = iSrc[i];
		}

		return srcByteSize;
	#endif

	#endif
}

uint64_t Uncompress(void * dst, void const * src, uint64_t srcByteSize, uint64_t dstMaxByteSize)
{
	#if K_USE_LZ4_COMPRESSION
	return LZ4_decompress_safe((char const *)src, (char *)dst, srcByteSize, dstMaxByteSize);
	#else
	// Dummy decompression test function
	#if 1
		uint32_t const * iSrc = (uint32_t const *)src;
		uint64_t * iDst = (uint64_t *)dst;

		const uint64_t n = srcByteSize / sizeof(uint32_t);
		for(uint64_t i = 0; i < n; ++i)
		{
			iDst[i] = (uint64_t)iSrc[i];
		}

		return srcByteSize * 2;
	#else
		uint64_t const * iSrc = (uint64_t const *)src;
		uint64_t * iDst = (uint64_t *)dst;

		const uint64_t n = srcByteSize / sizeof(uint64_t);
		for(uint64_t i = 0; i < n; ++i)
		{
			iDst[i] = iSrc[i];
		}

		return srcByteSize;
	#endif
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


// https://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#asynchronous-concurrent-execution
// https://www.docenti.unina.it/webdocenti-be/allegati/materiale-didattico/517066
/*
Two commands from different streams cannot run concurrently if one of the following operations is issued
in-between them by the host thread:
	- A page-locked host memory allocation,
	- A device memory allocation,
	- A device memory set,
	- A device <-> device memory copy,
	- A CUDA command to Stream 0 (including kernel launches and host <-> device memory copies that do
	  not specify any stream parameter)
	- A switch between the L1/shared memory configurations

	All GPU calls (memcpy, kernel execution, etc.) are placed into default stream unless otherwise specified
	
	Stream 0 is special:
		- Synchronous with all streams
		Meaning: Things done in Stream 0 cannot overlap with other streams
		- Streams with non-blocking flag are an exception:
			cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking)
*/

int PinnedHostMemoryTests()
{
	if(0)
	{
		const uint64_t byteSize = 8ull * 1024 * 1024 * 1024;
		unsigned char * pinnedHostMemory; 
		CUDA_CHECK_CALL(cudaHostAlloc(&pinnedHostMemory, byteSize, cudaHostAllocDefault));

		memset(pinnedHostMemory, 0xAB, byteSize);
		memset(pinnedHostMemory, 0xCD, byteSize);


		for(uint64_t i = 0; i < byteSize; ++i)
		{
			if(pinnedHostMemory[i] != 0xCD)
			{
				fprintf(stderr, "Incorrect value at %llu: 0x%x\n", i, (uint32_t)pinnedHostMemory[i]);
			}
		}
		CUDA_CHECK_CALL(cudaFreeHost(pinnedHostMemory));
	}

	if(0)
	{
		const uint64_t byteSize = 8ull * 1024 * 1024 * 1024;

		unsigned char * pinnedHostMemory; 
		CUDA_CHECK_CALL(cudaHostAlloc(&pinnedHostMemory, byteSize, cudaHostAllocDefault));

		unsigned char * hostMemory = (unsigned char *)malloc(byteSize);

		void * deviceMemory;
		CUDA_CHECK_CALL(cudaMalloc(&deviceMemory, byteSize));

		cudaStream_t stream;
		CUDA_CHECK_CALL(cudaStreamCreate(&stream));

		memset(pinnedHostMemory, 0xAB, byteSize);
		memset(pinnedHostMemory, 0xCD, byteSize);

		CUDA_CHECK_CALL(cudaMemcpyAsync(deviceMemory, pinnedHostMemory, byteSize, cudaMemcpyHostToDevice, stream));
		CUDA_CHECK_CALL(cudaStreamSynchronize(stream));

		CUDA_CHECK_CALL(cudaMemcpy(hostMemory, deviceMemory, byteSize, cudaMemcpyDeviceToHost));
		for(uint64_t i = 0; i < byteSize; ++i)
		{
			if(hostMemory[i] != 0xCD)
			{
				fprintf(stderr, "Incorrect value at %llu: 0x%x\n", i, (uint32_t)hostMemory[i]);
			}
		}

		CUDA_CHECK_CALL(cudaStreamDestroy(stream));
		CUDA_CHECK_CALL(cudaFree(deviceMemory));
		free(hostMemory);
		CUDA_CHECK_CALL(cudaFreeHost(pinnedHostMemory));
	}

	if(1)
	{
		const uint32_t numChunks = 4;
		const uint64_t chunkByteSize = 4ull * 1024ull * 1024ull; // 4MB
		const uint64_t byteSize = numChunks * chunkByteSize;

		unsigned char * deviceMemory; 
		CUDA_CHECK_CALL(cudaMalloc(&deviceMemory, byteSize));

		const uint32_t numStreams = 2;
		cudaStream_t memcpyStreams[numStreams];
		for(uint32_t i = 0; i < numStreams; ++i)
		{
			CUDA_CHECK_CALL(cudaStreamCreate(&memcpyStreams[i]));
		}

		unsigned char * chunkStagingData[numStreams];
		for(uint32_t i = 0; i < numStreams; ++i)
		{
			CUDA_CHECK_CALL(cudaHostAlloc(&chunkStagingData[i], chunkByteSize, cudaHostAllocDefault));
		}

		uint32_t currStreamIdx = 0;
		for(uint64_t i = 0; i < numChunks; ++i)
		{
			CUDA_CHECK_CALL(cudaStreamSynchronize(memcpyStreams[currStreamIdx])); // Ensure previous async copy is done

			unsigned char * src = chunkStagingData[currStreamIdx];
			unsigned char * dst = deviceMemory + i * chunkByteSize;

			memset(src, (currStreamIdx+1) * 0x11, chunkByteSize);
			memset(src, 0xAB, chunkByteSize);

			CUDA_CHECK_CALL(cudaMemcpyAsync(dst, src, chunkByteSize, cudaMemcpyHostToDevice, memcpyStreams[currStreamIdx])); // Not working

			currStreamIdx = (currStreamIdx + 1) % numStreams;
		}

		for(uint32_t i = 0; i < numStreams; ++i)
		{
			CUDA_CHECK_CALL(cudaStreamSynchronize(memcpyStreams[i])); // Ensure async copy is done
			CUDA_CHECK_CALL(cudaStreamDestroy(memcpyStreams[i]));;
		}

		for(uint32_t i = 0; i < numStreams; ++i)
		{
			CUDA_CHECK_CALL(cudaFreeHost(chunkStagingData[i]));
		}

		unsigned char * hostMemory = (unsigned char *)malloc(byteSize);
		CUDA_CHECK_CALL(cudaMemcpy(hostMemory, deviceMemory, byteSize, cudaMemcpyDeviceToHost));

		for(uint64_t i = 0; i < byteSize; ++i)
		{
			if(hostMemory[i] != 0xAB)
			{
				printf("Incorrect value %llu: %X\n", i, (uint32_t)hostMemory[i]);
			}
		}

		CUDA_CHECK_CALL(cudaFree(deviceMemory));
	}

	return 0;
}

// https://learn.microsoft.com/en-us/windows-hardware/drivers/display/device-paging-queues
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

	K_PUSH_PROFILING_MARKER(0xFFFFFFFF, "BinaryFormatTests");

	K_PUSH_PROFILING_MARKER(0xFFFF0000, "Init data");

	const uint64_t N = (1ull * 1024 * 1024 * 1024) / sizeof(uint64_t);
	const uint64_t srcDataByteSize = N * sizeof(uint64_t);
	uint64_t * srcData = (uint64_t *)malloc(srcDataByteSize);
	for(uint64_t i = 0; i < N; ++i)
		srcData[i] = i;

	K_POP_PROFILING_MARKER();

	K_PIX_BEGIN_EVENT(0xFFFF8800, "Allocate compression buffer");

	const uint64_t compressedBinaryDataByteSizeUpperBound = GetCompressedBinaryDataByteSizeUpperBound(srcDataByteSize);
	void * compressionBuffer = malloc(compressedBinaryDataByteSizeUpperBound);
	
	K_POP_PROFILING_MARKER();

	K_PUSH_PROFILING_MARKER(0xFFFFFF00, "Compress data");

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
	
	K_POP_PROFILING_MARKER();

	const uint64_t uncompressedBinaryDataByteSize = GetUncompressedBinaryDataByteSize(compressionBuffer);
	assert(uncompressedBinaryDataByteSize == srcDataByteSize);

	void * decompressionBuffer = malloc(uncompressedBinaryDataByteSize);

	const auto ValidateDecompressionData = [srcData, srcDataByteSize](char const * label, void * decompressionBuffer)
	{
		uint64_t const * uncompressedData = (uint64_t const *)decompressionBuffer;
		int cmp = memcmp(srcData, uncompressedData, srcDataByteSize);
		if(cmp != 0)
		{
			fprintf(stderr, "Error while decompressing: uncompressedData != srcData [%s]\n", label);
		}
	};

	// Single-thread - decompress into local buffer + upload to CUDA buffer
	if(1)
	{
		K_PUSH_PROFILING_MARKER(0xFFFFFF00, "Device memory - alloc");
		void * deviceMemory = nullptr;
		CUDA_CHECK_CALL(cudaMalloc(&deviceMemory, srcDataByteSize));
		K_POP_PROFILING_MARKER();

		uint64_t uncompressedBinaryDataByteSize2 = 0;

		double memsetHostMemoryDuration_ms;
		double decompressionDuration_ms;
		double memcpyHostToDevice_ms;
		MeasureAndPrintFuncTime("Single-thread decompression", [&]()
		{
			memsetHostMemoryDuration_ms = MeasureFuncTime([&]()
			{
				K_PUSH_PROFILING_MARKER(0xFFFFFF00, "Memset host data");
				memset(decompressionBuffer, 0xFF, uncompressedBinaryDataByteSize);
				K_POP_PROFILING_MARKER();
			});

			decompressionDuration_ms = MeasureFuncTime([&]()
			{
				K_PUSH_PROFILING_MARKER(0xFF00FF00, "Uncompress data - standard");
				uncompressedBinaryDataByteSize2 = UncompressedBinaryData(decompressionBuffer, compressionBuffer, compressedBinaryDataByteSize, uncompressedBinaryDataByteSize);
				K_POP_PROFILING_MARKER();
			});
			
			memcpyHostToDevice_ms = MeasureFuncTime([&]()
			{
				K_PUSH_PROFILING_MARKER(0xFF00FFFF, "Memcpy host to device");
				CUDA_CHECK_CALL(cudaMemcpy(deviceMemory, decompressionBuffer, srcDataByteSize, cudaMemcpyHostToDevice));
				K_POP_PROFILING_MARKER();
			});
		});
		printf("  Memset host memory: %fms\n", memsetHostMemoryDuration_ms);
		printf("  Decompression: %fms\n", decompressionDuration_ms);
		printf("  Memcpy host -> device: %fms (bandwidth: %s/s)\n", memcpyHostToDevice_ms, HumanReadableByteSizeStr(srcDataByteSize / (memcpyHostToDevice_ms / 1000.0)).c_str());

		assert(uncompressedBinaryDataByteSize2 == uncompressedBinaryDataByteSize);
		assert(uncompressedBinaryDataByteSize2 == srcDataByteSize);

		ValidateDecompressionData("Single-thread", decompressionBuffer);

		K_PUSH_PROFILING_MARKER(0xFF0000FF, "Device memory - free");
		CUDA_CHECK_CALL(cudaFree(deviceMemory));
		K_POP_PROFILING_MARKER();
	}

	// Parallel - decompress into local buffer + upload to CUDA buffer
	if(0)
	{
		K_PUSH_PROFILING_MARKER(0xFF0000FF, "Uncompress data - parallel - CUDA");
		memset(decompressionBuffer, 0xFF, uncompressedBinaryDataByteSize);

		struct DataChunkToProcess
		{
			unsigned char * dst;
			unsigned char const * src;
			uint32_t srcByteSize;
			uint32_t dstByteSize;
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

		Barrier barrier(numThreads + 1);
		Barrier & initializationBarrier = barrier;
		Barrier & processingDoneBarrier = barrier;

		for(int i = 0; i < numThreads; ++i)
		{
			static constexpr uint32_t kNumBufferedChunks = 12;

			#define USE_WRITE_COMBINED_MEMORY 0 // Note: **Really** slow because the LZ4 decompression routine read from the output buffer while decompressing

			#if USE_WRITE_COMBINED_MEMORY
			const unsigned int hostAllocFlags = cudaHostAllocWriteCombined;
			#else
			const unsigned int hostAllocFlags = cudaHostAllocDefault;
			#endif

			const uint64_t tmpChunkStorageByteSize = kNumBufferedChunks * BinDataHeader::MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE;
			unsigned char * tmpChunkPinnedHostMemoryStorage;

			// https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html#asynchronous-and-overlapping-transfers-with-computation
			CUDA_CHECK_CALL(cudaHostAlloc(&tmpChunkPinnedHostMemoryStorage, tmpChunkStorageByteSize, hostAllocFlags));
			CUDA_CHECK_CALL(cudaMemset(tmpChunkPinnedHostMemoryStorage, 0x00, tmpChunkStorageByteSize)); // Warm up pinned memory ("MakeResident" + "CommitVirtualAddressRange" Paging Queue Packets on the "UMD Paging queue" in NVIDIA NSight Systems)

			workers[i] = std::thread([&, tmpChunkPinnedHostMemoryStorage]()
			{
				unsigned char * tmpChunkPinnedHostMemoryBuffer[kNumBufferedChunks];
				for(uint32_t j = 0; j < kNumBufferedChunks; ++j)
				{
					tmpChunkPinnedHostMemoryBuffer[j] = tmpChunkPinnedHostMemoryStorage + j * BinDataHeader::MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE;
				}

				cudaStream_t memcpyStreams[kNumBufferedChunks];
				for(int j = 0; j < kNumBufferedChunks; ++j)
				{
					CUDA_CHECK_CALL(cudaStreamCreate(&memcpyStreams[j]));
					//CUDA_CHECK_CALL(cudaEventCreateWithFlags(&memcpyEvents[j], cudaEventBlockingSync|cudaEventDisableTiming));
				}

				initializationBarrier.wait(); // Wait for initialization of all workers and main thread

				uint32_t currBuffer = 0;
				while(true)
				{
					K_PUSH_PROFILING_MARKER(0xFFFF00FF, "Uncompress chunk - wait");
					std::unique_lock<std::mutex> lock(queueMutex);
					queueCV.wait(lock, [&]{ return (/*bReady && */!dataChunksToProcess.empty()) || bKillThreads; });

					if(bKillThreads)
						break;
					K_POP_PROFILING_MARKER();

					K_PUSH_PROFILING_MARKER(0xFF0077FF, "Uncompress chunk - pop item");
					DataChunkToProcess chunk = dataChunksToProcess.front();
					dataChunksToProcess.pop();
					lock.unlock();
					K_POP_PROFILING_MARKER();

					K_PUSH_PROFILING_MARKER(0xFF0000FF, "cudaStreamSynchronize");
					CUDA_CHECK_CALL(cudaStreamSynchronize(memcpyStreams[currBuffer])); // Ensure previous async copy is done
					K_POP_PROFILING_MARKER();

					K_PUSH_PROFILING_MARKER(0xFF00FFFF, "Uncompress chunk - parallel - CUDA");
					assert(chunk.dstByteSize <= BinDataHeader::MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE);
					Uncompress(tmpChunkPinnedHostMemoryBuffer[currBuffer], chunk.src, chunk.srcByteSize, chunk.dstByteSize);
					K_POP_PROFILING_MARKER();

					K_PUSH_PROFILING_MARKER(0xFF0000FF, "cudaMemcpyAsync");
					CUDA_CHECK_CALL(cudaMemcpyAsync(chunk.dst, tmpChunkPinnedHostMemoryBuffer[currBuffer], chunk.dstByteSize, cudaMemcpyHostToDevice, memcpyStreams[currBuffer]));
					//CUDA_CHECK_CALL(cudaMemcpy(chunk.dst, tmpChunkPinnedHostMemoryBuffer[currBuffer], chunk.dstByteSize, cudaMemcpyHostToDevice));
					K_POP_PROFILING_MARKER();

					currBuffer = (currBuffer + 1) % kNumBufferedChunks;
				}

				for(int j = 0; j < kNumBufferedChunks; ++j)
				{
					CUDA_CHECK_CALL(cudaStreamSynchronize(memcpyStreams[j])); // Ensure async copy is done
				}

				processingDoneBarrier.wait(); // Wait for workers exiting their processing loop

				// Clean up
				for(int j = 0; j < kNumBufferedChunks; ++j)
				{
					CUDA_CHECK_CALL(cudaStreamDestroy(memcpyStreams[j]));;
				}
				CUDA_CHECK_CALL(cudaFreeHost(tmpChunkPinnedHostMemoryStorage));
			});
		}

		void * decompressionBufferDevice;
		CUDA_CHECK_CALL(cudaMalloc(&decompressionBufferDevice, uncompressedBinaryDataByteSize));

		K_PUSH_PROFILING_MARKER(0xFF00FFFF, "Parallel decompression and CUDA buffer upload");
		start = std::chrono::high_resolution_clock::now();

		const uint64_t numChunks = GetCompressedBinaryDataNumChunks(compressionBuffer);
		//dataChunksToProcess.reserve(numChunks);

		struct Userdata
		{
			void * pBaseDstChunkDataDevice;
			std::queue<DataChunkToProcess> & dataChunksToProcess;
		};

		Userdata userdata = { decompressionBufferDevice, dataChunksToProcess };
		const auto chunkDecompressionFunc = [](unsigned char const * pSrcChunkData, uint32_t srcChunkByteSize, uint64_t dstByteOffset, uint32_t dstChunkByteSize, void * pUserdata)
		{
			Userdata const & userdata = *(Userdata const *)pUserdata;
			unsigned char * pDstChunkData = (unsigned char*)userdata.pBaseDstChunkDataDevice + dstByteOffset;
			userdata.dataChunksToProcess.push({pDstChunkData, pSrcChunkData, srcChunkByteSize, dstChunkByteSize});
		};

		K_PUSH_PROFILING_MARKER(0xFF00FFFF, "Prepare chunks - parallel - CUDA");
		const uint64_t uncompressedBinaryDataByteSize2 = UncompressedBinaryData(compressionBuffer, compressedBinaryDataByteSize, chunkDecompressionFunc, &userdata);
		K_POP_PROFILING_MARKER();

		initializationBarrier.wait(); // Wait for workers initialization (TODO: add thread-safe queue to enqueue and process in parallel)

		#if 0
		while(!dataChunksToProcess.empty())
		{
			DataChunkToProcess const & chunk = dataChunksToProcess.front();
			Uncompress(chunk.dst, chunk.src, chunk.srcByteSize, chunk.dstMaxByteSize);
			dataChunksToProcess.pop();
		}
		#else
		// Notify workers to start working
		//{
		//	std::unique_lock<std::mutex> lock(queueMutex);
		//	bReady = true;
		//}
		queueCV.notify_all();
		#endif

		// Waits for workers to finish
		// TODO: find a better way to check that workers are done
		K_PUSH_PROFILING_MARKER(0xFFFF00FF, "Main thread - wait for workers"); // TODO: make main thread participate as well?
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

			std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}
		K_POP_PROFILING_MARKER();

		processingDoneBarrier.wait(); // Wait for workers exiting their processing loop

		K_POP_PROFILING_MARKER();

		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		duration_ms = duration.count() / 1000.0;
		printf("Decompression time and GPU upload time: %fms - [parallel - CUDA]\n", duration_ms);

		for(std::thread & w : workers)
			w.join();

		assert(uncompressedBinaryDataByteSize2 == uncompressedBinaryDataByteSize);
		assert(uncompressedBinaryDataByteSize2 == srcDataByteSize);

		CUDA_CHECK_CALL(cudaMemcpy(decompressionBuffer, decompressionBufferDevice, uncompressedBinaryDataByteSize, cudaMemcpyDeviceToHost));

		uint64_t const * uncompressedData = (uint64_t const *)decompressionBuffer;

		unsigned char const * uncompressedDataU8 = (unsigned char const *)decompressionBuffer;

		ValidateDecompressionData("Multi-thread", decompressionBuffer);

		CUDA_CHECK_CALL(cudaFree(decompressionBufferDevice));

	
	if(1)
	{
		K_PUSH_PROFILING_MARKER(0xFF0000FF, "Uncompress data - parallel - CUDA");
		memset(decompressionBuffer, 0xFF, uncompressedBinaryDataByteSize);

		struct DataChunkToProcess
		{
			unsigned char * dst;
			unsigned char const * src;
			uint32_t srcByteSize;
			uint32_t dstByteSize;
		};

		// TODO: use SPMC queue?
		std::queue<DataChunkToProcess> dataChunksToProcess;

		struct DecompressionWorkspace
		{
			unsigned char * tmpChunkPinnedHostMemoryBuffer;
			size_t tmpChunkPinnedHostMemoryByteSize;
			cudaStream_t memcpyStream;
			std::mutex tmpDecompressionBufferMutex;
			std::condition_variable tmpDecompressionBufferCV;
			bool bDecompressing = false;
		};

		std::mutex queueMutex;
		std::condition_variable queueCV;
		bool bReady = false;
		bool bKillThreads = false;

		const int numThreads = std::thread::hardware_concurrency() - 2;
		printf("Num threads: %d\n", numThreads);
		std::vector<std::thread> workers(numThreads);

		#if 0
		static constexpr uint32_t kNumBufferedChunks = 12;
		const uint64_t numDecompressionWorkspace = numThreads * kNumBufferedChunks;
		#else
		const uint64_t maxPinnedHostMemoryByteSize = 2ull * 1024ull * 1024ull * 1024ull;
		const uint64_t numDecompressionWorkspace = maxPinnedHostMemoryByteSize / BinDataHeader::MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE;
		// Note: would also need to limit the number of threads if numDecompressionWorkspace < thread_pool_size * some_multiplier
		#endif

		std::vector<DecompressionWorkspace> decompressionWS(numDecompressionWorkspace);

		#define USE_WRITE_COMBINED_MEMORY 0 // Note: **Really** slow because the LZ4 decompression routine read from the output buffer while decompressing

		#if USE_WRITE_COMBINED_MEMORY
		const unsigned int hostAllocFlags = cudaHostAllocWriteCombined;
		#else
		const unsigned int hostAllocFlags = cudaHostAllocDefault;
		#endif

		const uint64_t decompressionWSTotalByteSize = decompressionWS.size() * BinDataHeader::MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE;
		// https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html#asynchronous-and-overlapping-transfers-with-computation
		unsigned char * decompressionWSPinnedHostMemoryStorage;
		CUDA_CHECK_CALL(cudaHostAlloc(&decompressionWSPinnedHostMemoryStorage, decompressionWSTotalByteSize, hostAllocFlags));
		//CUDA_CHECK_CALL(cudaMemset(decompressionWSPinnedHostMemoryStorage, 0x00, decompressionWSTotalByteSize)); // Warm up pinned memory ("MakeResident" + "CommitVirtualAddressRange" Paging Queue Packets on the "UMD Paging queue" in NVIDIA NSight Systems)
		
		for(size_t i = 0; i < decompressionWS.size(); ++i)
		{
			decompressionWS[i].tmpChunkPinnedHostMemoryBuffer = decompressionWSPinnedHostMemoryStorage + i * BinDataHeader::MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE;
			decompressionWS[i].tmpChunkPinnedHostMemoryByteSize = BinDataHeader::MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE;
			CUDA_CHECK_CALL(cudaStreamCreate(&decompressionWS[i].memcpyStream));
		}

		Barrier barrier(numThreads + 1);
		Barrier & initializationBarrier = barrier;
		Barrier & processingDoneBarrier = barrier;

		std::atomic<uint32_t> decompressionWSIndex = 0;

		for(int i = 0; i < numThreads; ++i)
		{
			workers[i] = std::thread([&]()
			{
				initializationBarrier.wait(); // Wait for initialization of all workers and main thread

				while(true)
				{
					K_PUSH_PROFILING_MARKER(0xFFFF00FF, "Uncompress chunk - wait");
					std::unique_lock<std::mutex> lock(queueMutex);
					queueCV.wait(lock, [&]{ return (/*bReady && */!dataChunksToProcess.empty()) || bKillThreads; });
					K_POP_PROFILING_MARKER();

					if(bKillThreads)
						break;

					K_PUSH_PROFILING_MARKER(0xFF0077FF, "Uncompress chunk - pop item");
					DataChunkToProcess chunk = dataChunksToProcess.front();
					dataChunksToProcess.pop();
					lock.unlock();
					K_POP_PROFILING_MARKER();

					uint32_t index = decompressionWSIndex++;
					index = index % decompressionWS.size();

					K_PUSH_PROFILING_MARKER(0xFF885522, "DecompressionWorkspace - wait");
					DecompressionWorkspace & decompressionWorkspace = decompressionWS[index];
					std::unique_lock<std::mutex> decompressionWSLock(decompressionWorkspace.tmpDecompressionBufferMutex);
					decompressionWorkspace.tmpDecompressionBufferCV.wait(decompressionWSLock, [&]{ return !decompressionWorkspace.bDecompressing; });
					K_POP_PROFILING_MARKER();

					decompressionWorkspace.bDecompressing = true;

					K_PUSH_PROFILING_MARKER(0xFF0000FF, "cudaStreamSynchronize");
					CUDA_CHECK_CALL(cudaStreamSynchronize(decompressionWorkspace.memcpyStream)); // Ensure previous async copy is done
					K_POP_PROFILING_MARKER();

					K_PUSH_PROFILING_MARKER(0xFF00FFFF, "Uncompress chunk - parallel - CUDA");
					assert(chunk.dstByteSize <= BinDataHeader::MAX_CHUNKS_UNCOMPRESSED_BYTE_SIZE);
					Uncompress(decompressionWorkspace.tmpChunkPinnedHostMemoryBuffer, chunk.src, chunk.srcByteSize, chunk.dstByteSize);
					K_POP_PROFILING_MARKER();

					K_PUSH_PROFILING_MARKER(0xFF0000FF, "cudaMemcpyAsync");
					CUDA_CHECK_CALL(cudaMemcpyAsync(chunk.dst, decompressionWorkspace.tmpChunkPinnedHostMemoryBuffer, chunk.dstByteSize, cudaMemcpyHostToDevice, decompressionWorkspace.memcpyStream));
					//CUDA_CHECK_CALL(cudaMemcpy(chunk.dst, decompressionWorkspace.tmpChunkPinnedHostMemoryBuffer, chunk.dstByteSize, cudaMemcpyHostToDevice));
					K_POP_PROFILING_MARKER();

					K_PUSH_PROFILING_MARKER(0xFFFFFF00, "DecompressionWorkspace - unlock and notify");
					decompressionWorkspace.bDecompressing = false;
					decompressionWSLock.unlock();
					decompressionWorkspace.tmpDecompressionBufferCV.notify_one();
					K_POP_PROFILING_MARKER();
				}

				processingDoneBarrier.wait(); // Wait for workers exiting their processing loop
			});
		}

		void * decompressionBufferDevice;
		CUDA_CHECK_CALL(cudaMalloc(&decompressionBufferDevice, uncompressedBinaryDataByteSize));

		K_PUSH_PROFILING_MARKER(0xFF00FFFF, "Parallel decompression and CUDA buffer upload");
		start = std::chrono::high_resolution_clock::now();

		const uint64_t numChunks = GetCompressedBinaryDataNumChunks(compressionBuffer);
		//dataChunksToProcess.reserve(numChunks);

		struct Userdata
		{
			void * pBaseDstChunkDataDevice;
			std::queue<DataChunkToProcess> & dataChunksToProcess;
		};

		Userdata userdata = { decompressionBufferDevice, dataChunksToProcess };
		const auto chunkDecompressionFunc = [](unsigned char const * pSrcChunkData, uint32_t srcChunkByteSize, uint64_t dstByteOffset, uint32_t dstChunkByteSize, void * pUserdata)
		{
			Userdata const & userdata = *(Userdata const *)pUserdata;
			unsigned char * pDstChunkData = (unsigned char*)userdata.pBaseDstChunkDataDevice + dstByteOffset;
			userdata.dataChunksToProcess.push({pDstChunkData, pSrcChunkData, srcChunkByteSize, dstChunkByteSize});
		};

		K_PUSH_PROFILING_MARKER(0xFF00FFFF, "Prepare chunks - parallel - CUDA");
		const uint64_t uncompressedBinaryDataByteSize2 = UncompressedBinaryData(compressionBuffer, compressedBinaryDataByteSize, chunkDecompressionFunc, &userdata);
		K_POP_PROFILING_MARKER();

		initializationBarrier.wait(); // Wait for workers initialization (TODO: add thread-safe queue to enqueue and process in parallel)

		#if 0
		while(!dataChunksToProcess.empty())
		{
			DataChunkToProcess const & chunk = dataChunksToProcess.front();
			Uncompress(chunk.dst, chunk.src, chunk.srcByteSize, chunk.dstMaxByteSize);
			dataChunksToProcess.pop();
		}
		#else
		// Notify workers to start working
		//{
		//	std::unique_lock<std::mutex> lock(queueMutex);
		//	bReady = true;
		//}
		queueCV.notify_all();
		#endif

		// Waits for workers to finish
		// TODO: find a better way to check that workers are done
		K_PUSH_PROFILING_MARKER(0xFFFF00FF, "Main thread - wait for workers"); // TODO: make main thread participate as well?
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

			std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}
		K_POP_PROFILING_MARKER();

		K_PUSH_PROFILING_MARKER(0xFF55FF88, "cudaDeviceSynchronize"); // TODO: make main thread participate as well?
		CUDA_CHECK_CALL(cudaDeviceSynchronize()); // Wait for all copies to finish (TODO: try to call cudaStreamSynchronize in parallel by worker threads?)
		K_POP_PROFILING_MARKER();

		processingDoneBarrier.wait(); // Wait for workers exiting their processing loop

		K_POP_PROFILING_MARKER();

		stop = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		duration_ms = duration.count() / 1000.0;
		printf("Decompression time and GPU upload time: %fms - [parallel - CUDA]\n", duration_ms);

		K_PUSH_PROFILING_MARKER(0xFFFF7766, "Workers join");
		for(std::thread & w : workers)
			w.join();
		K_POP_PROFILING_MARKER();

		K_PUSH_PROFILING_MARKER(0xFFFFFFFF, "Validate data");
		assert(uncompressedBinaryDataByteSize2 == uncompressedBinaryDataByteSize);
		assert(uncompressedBinaryDataByteSize2 == srcDataByteSize);

		CUDA_CHECK_CALL(cudaMemcpy(decompressionBuffer, decompressionBufferDevice, uncompressedBinaryDataByteSize, cudaMemcpyDeviceToHost));

		uint64_t const * uncompressedData = (uint64_t const *)decompressionBuffer;

		unsigned char const * uncompressedDataU8 = (unsigned char const *)decompressionBuffer;

		ValidateDecompressionData("Multi-thread", decompressionBuffer);
		K_POP_PROFILING_MARKER();

		K_PUSH_PROFILING_MARKER(0xFFFF5522, "Clean up");
		CUDA_CHECK_CALL(cudaFree(decompressionBufferDevice));

		for(size_t i = 0; i < decompressionWS.size(); ++i)
		{
			CUDA_CHECK_CALL(cudaStreamDestroy(decompressionWS[i].memcpyStream));
		}
		CUDA_CHECK_CALL(cudaFreeHost(decompressionWSPinnedHostMemoryStorage));
		K_POP_PROFILING_MARKER();

		K_POP_PROFILING_MARKER();
	}

	/*
	   // https://stackoverflow.com/questions/49480334/where-is-pinned-memory-allocated-using-cudahostalloc

		"Page-Locked Host Memory" for CUDA (and other DMA-capable external hardware like PCI-express cards) is allocated in physical memory of the Host computer. The allocation is marked as not-swappable (not-pageable) and not-movable (locked, pinned). This is similar to the action of mlock syscall "lock part or all of the calling process's virtual address space into RAM, preventing that memory from being paged to the swap area."

		This allocation can be accessed by kernel virtual address space (as kernel has full view of the physical memory) and this allocation is also added to the user process virtual address space to allow process access it.

		When you does ordinary malloc, actual physical memory allocation may (and will) be postponed to the first (write) access to the pages. With mlocked/pinned memory all physical pages are allocated inside locking or pinning calls (like MAP_POPULATE in mmap: "Populate (prefault) page tables for a mapping"), and physical addresses of pages will not change (no swapping, no moving, no compacting...).

		CUDA docs: http://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html#group__CUDART__MEMORY_1gb65da58f444e7230d3322b6126bb4902

			__host__ cudaError_t cudaHostAlloc ( void** pHost, size_t size, unsigned int  flags )

			Allocates page-locked memory on the host. ...

			Allocates size bytes of host memory that is page-locked and accessible to the device. The driver tracks the virtual memory ranges allocated with this function and automatically accelerates calls to functions such as cudaMemcpy(). Since the memory can be accessed directly by the device, it can be read or written with much higher bandwidth than pageable memory obtained with functions such as malloc(). Allocating excessive amounts of pinned memory may degrade system performance, since it reduces the amount of memory available to the system for paging. As a result, this function is best used sparingly to allocate staging areas for data exchange between host and device.

			...

			Memory allocated by this function must be freed with cudaFreeHost().

		Pinned and not-pinned memory compared: https://www.cs.virginia.edu/~mwb7w/cuda_support/pinned_tradeoff.html "Choosing Between Pinned and Non-Pinned Memory"

			Pinned memory is memory allocated using the cudaMallocHost function, which prevents the memory from being swapped out and provides improved transfer speeds. Non-pinned memory is memory allocated using the malloc function. As described in Memory Management Overhead and Memory Transfer Overhead, pinned memory is much more expensive to allocate and deallocate but provides higher transfer throughput for large memory transfers.

		CUDA forums post with advises from txbob moderator: https://devtalk.nvidia.com/default/topic/899020/does-cudamemcpyasync-require-pinned-memory-/ "Does cudaMemcpyAsync require pinned memory?"

			If you want truly asynchronous behavior (e.g. overlap of copy and compute) then the memory must be pinned. If it is not pinned, there won't be any runtime errors, but the copy will not be asynchronous - it will be performed like an ordinary cudaMemcpy.

			The usable size may vary by system and OS. Pinning 4GB of memory on a 64GB system on Linux should not have a significant effect on CPU performance, after the pinning operation is complete. Attempting to pin 60GB on the other hand might cause significant system responsiveness issues.
	*/

	// Memcpy - from default pinned host memory
	if(0)
	{
		K_PUSH_PROFILING_MARKER(0xFFFFFF00, "Pinned host memory - alloc");
		void * pinnedHostMemory = nullptr;
		CUDA_CHECK_CALL(cudaHostAlloc(&pinnedHostMemory, srcDataByteSize, cudaHostAllocDefault));
		K_POP_PROFILING_MARKER();
		
		K_PUSH_PROFILING_MARKER(0xFFFFFF00, "Device memory - alloc");
		void * deviceMemory = nullptr;
		CUDA_CHECK_CALL(cudaMalloc(&deviceMemory, srcDataByteSize));
		K_POP_PROFILING_MARKER();

		const auto start = std::chrono::high_resolution_clock::now();
		memcpy(pinnedHostMemory, srcData, srcDataByteSize);
		CUDA_CHECK_CALL(cudaMemcpy(deviceMemory, pinnedHostMemory, srcDataByteSize, cudaMemcpyHostToDevice));
		const auto stop = std::chrono::high_resolution_clock::now();
		const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
		const double duration_ms = duration.count() / 1000.0;
		printf("Decompression time: %fms - [memcpy pinned host memory]\n", duration_ms);

		K_PUSH_PROFILING_MARKER(0xFF0000FF, "Pinned host memory - free");
		CUDA_CHECK_CALL(cudaFreeHost(pinnedHostMemory));
		K_POP_PROFILING_MARKER();
		
		K_PUSH_PROFILING_MARKER(0xFF0000FF, "Device memory - free");
		CUDA_CHECK_CALL(cudaFree(deviceMemory));
		K_POP_PROFILING_MARKER();
	}
	
	// Memcpy - from write-combined pinned host memory
	if(1)
	{
		K_PUSH_PROFILING_MARKER(0xFFFFFF00, "Pinned host memory - alloc");
		void * pinnedHostMemory = nullptr;
		CUDA_CHECK_CALL(cudaHostAlloc(&pinnedHostMemory, srcDataByteSize, cudaHostAllocWriteCombined));
		K_POP_PROFILING_MARKER();
		
		K_PUSH_PROFILING_MARKER(0xFFFFFF00, "Device memory - alloc");
		void * deviceMemory = nullptr;
		CUDA_CHECK_CALL(cudaMalloc(&deviceMemory, srcDataByteSize));
		K_POP_PROFILING_MARKER();

		double memcpyHostToPinnedHostMemoryDuration_ms;
		double memcpyPinnedHostMemoryToDeviceDuration_ms;
		MeasureAndPrintFuncTime("Memcpy - from write-combined pinned host memory", [&]()
		{
			memcpyHostToPinnedHostMemoryDuration_ms = MeasureFuncTime([&]()
			{
				K_PUSH_PROFILING_MARKER(0xFF0000FF, "Memcpy: Host -> Pinned host WC memory");
				memcpy(pinnedHostMemory, srcData, srcDataByteSize);
				K_POP_PROFILING_MARKER();
			});

			memcpyPinnedHostMemoryToDeviceDuration_ms = MeasureFuncTime([&]()
			{
				K_PUSH_PROFILING_MARKER(0xFF0000FF, "Memcpy: Pinned host WC memory -> Device");
				CUDA_CHECK_CALL(cudaMemcpy(deviceMemory, pinnedHostMemory, srcDataByteSize, cudaMemcpyHostToDevice));
				K_POP_PROFILING_MARKER();
			});
		});
		printf("  Memcpy: host -> pinned host memory: %fms (%s/s)\n", memcpyHostToPinnedHostMemoryDuration_ms, HumanReadableByteSizeStr(srcDataByteSize / (memcpyHostToPinnedHostMemoryDuration_ms / 1000.0)).c_str());
		printf("  Memcpy: pinned host memory -> device: %fms (%s/s)\n", memcpyPinnedHostMemoryToDeviceDuration_ms, HumanReadableByteSizeStr(srcDataByteSize / (memcpyPinnedHostMemoryToDeviceDuration_ms / 1000.0)).c_str());

		K_PUSH_PROFILING_MARKER(0xFF0000FF, "Pinned host memory - free");
		CUDA_CHECK_CALL(cudaFreeHost(pinnedHostMemory));
		K_POP_PROFILING_MARKER();
		
		K_PUSH_PROFILING_MARKER(0xFF0000FF, "Device memory - free");
		CUDA_CHECK_CALL(cudaFree(deviceMemory));
		K_POP_PROFILING_MARKER();
	}

	free(srcData);
	free(compressionBuffer);
	free(decompressionBuffer);

	K_POP_PROFILING_MARKER();

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
	CUDA_DRIVER_CHECK_CALL(cuInit(0));
	CUDA_DRIVER_CHECK_CALL(cuDeviceGet(&cuDevice, 0));
	CUDA_DRIVER_CHECK_CALL(cuCtxCreate(&context, 0, cuDevice));

	std::mutex mutex;
	std::condition_variable cv;
	bool bReady = false;

	const int numThreads = std::thread::hardware_concurrency() - 2;
	std::vector<std::thread> workers(numThreads);
	for(int i = 0; i < numThreads; ++i)
	{
		workers[i] = std::thread([&, i]()
		{
			CUDA_DRIVER_CHECK_CALL(cuCtxSetCurrent(context));

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

	CUDA_DRIVER_CHECK_CALL(cuCtxDestroy(context));
	#endif

	return 0;
}

#pragma endregion