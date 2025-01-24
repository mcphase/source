#include "ic1ion_mod.hpp"
#if defined(__linux__) || defined(__APPLE__)
extern "C"
{
ic1ion_mod *allocator()
{
return new ic1ion_mod();
}
void deleter(ic1ion_mod *ptr)
{
delete ptr;
}
}
#endif
#ifdef WIN32
extern "C"
{
__declspec (dllexport) ic1ion_mod *allocator()
{
return new ic1ion_mod();
}
__declspec (dllexport) void deleter(ic1ion_mod *ptr)
{
delete ptr;
}
}
#endif


void ic1ion_mod::Icalc()
{
std::cout << "Greetings from ic1ion_mod !" << std::endl;
}

