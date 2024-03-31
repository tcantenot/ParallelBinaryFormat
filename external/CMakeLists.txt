project(external C CXX)

# LZ4 compressionlib/
add_library(lz4 EXCLUDE_FROM_ALL
	lz4/lib/lz4.h
	lz4/lib/lz4.c
)
target_include_directories(lz4 PUBLIC lz4/lib)
add_library(ext::lz4 ALIAS lz4)
