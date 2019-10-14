#include "utils.h"

#include <time.h>
#include <iostream>
#include <string>
#include <cstdio>
// Utils for measuring time...

//Returns the size (in bytes) of a file...
size_t getFileSize(FILE *fp) {	
    if (fp) {
        #ifndef _WIN32	
            fseek(fp,0,SEEK_END);
            size_t size = ftell(fp);
            rewind(fp);
        #else
            _fseeki64(fp, (__int64)0, SEEK_END);
        	size_t size = _ftelli64(fp);
        	rewind(fp);
        #endif
        return size;
    }
    return 0;
}

//reads a file and returns true iff everything went ok...
bool readFileContent(char *buffer, size_t sizeRead, FILE *fp)  {
    size_t sizeActuallyRead = fread(buffer,1,sizeRead,fp);
    return sizeActuallyRead!=0;
}


