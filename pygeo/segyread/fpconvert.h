#include <Python.h>
#include <sys/types.h>
#include <netinet/in.h>

void ibm2ieee (void *to, const void *from, Py_ssize_t len);
void ieee2ibm (void *to, const void *from, Py_ssize_t len);
