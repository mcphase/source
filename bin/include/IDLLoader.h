#pragma once
#include <memory>
#include <string>
namespace dlloader
{
template <class T>
class IDLLoader
{
public:
IDLLoader() = default;
virtual ~IDLLoader() = default;
/*
** Load the library and map it in memory.
*/
virtual void DLOpenLib() = 0;
/*
** Return a shared pointer on an instance of class loaded through
** a dynamic library.
*/
virtual std::shared_ptr<T>  DLGetInstance() = 0;
/*
** Unload the library.
*/
virtual void DLCloseLib() = 0;
};
}
