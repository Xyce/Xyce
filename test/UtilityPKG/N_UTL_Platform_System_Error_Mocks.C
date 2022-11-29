#include <Xyce_config.h>

#if defined(HAVE_WINDOWS_H)
#include <Windows.h>
#include <Sysinfoapi.h>

#elif defined(HAVE_UNISTD_H) && defined(HAVE_GETHOSTNAME)
#include <unistd.h>
#include <errno.h>
#endif

// hostname system calls
#if defined(HAVE_WINDOWS_H)
// TO DO - Mock implementation of GetComputerNameEx
// BOOL GetComputerNameEx(COMPUTER_NAME_FORMAT NameType, LPSTR lpBuffer, LPDWORD nSize);

#elif defined(HAVE_GETHOSTNAME)
int gethostname(char *__name, size_t __length) throw()
{
  errno = EFAULT;
  return -1;
}

#endif