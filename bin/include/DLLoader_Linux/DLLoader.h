#pragma once           //DLLoader.h for Linux
#include <iostream>
#include <dlfcn.h>
#include "IDLLoader.h"

namespace dlloader
{

template <typename T> class DLLoader : public IDLLoader<T>
{
private:

void *_handle;
std::string _pathToLib;
std::string _allocClassSymbol;
std::string _deleteClassSymbol;

public: 
DLLoader(std::string const &pathToLib,
         std::string const &allocClassSymbol = "allocator",
         std::string const &deleteClassSymbol = "deleter") : _handle(nullptr),
                    _pathToLib(pathToLib),_allocClassSymbol(allocClassSymbol),
                    _deleteClassSymbol(deleteClassSymbol){}

DLLoader() = default;
~DLLoader() = default;

void DLOpenLib() override
{
if (!(_handle = dlopen(_pathToLib.c_str(), RTLD_NOW | RTLD_LAZY))) {std::cerr << dlerror() << std::endl;exit(EXIT_FAILURE); }
}

void DLCloseLib() override
{
if (dlclose(_handle) != 0) {std::cerr << dlerror() << std::endl;}
}

std::shared_ptr <T> DLGetInstance() override{

//using allocClass = T *(*)();
//using deleteClass = void (*)(T *);

//void(*allocFunc)()  = reinterpret_cast<void(*)()>(dlsym(_handle, _allocClassSymbol.c_str()));
//void(*deleteFunc)() = reinterpret_cast<void(*)()>(dlsym(_handle, _deleteClassSymbol.c_str()));
// auto_ptr
// unique_ptr
//weak_ptr
// nullptr_t
// shared_ptr
auto allocFunc = reinterpret_cast<T *(*)()>(dlsym(_handle, _allocClassSymbol.c_str()));
auto deleteFunc = reinterpret_cast<void (*)(T *)>(dlsym(_handle, _deleteClassSymbol.c_str()));
if (!allocFunc || !deleteFunc) { std::cerr << dlerror() << std::endl;DLCloseLib();exit(EXIT_FAILURE);}

//return std::shared_ptr<T>(allocFunc(),[deleteFunc](T *p){ deleteFunc(p); });
return std::shared_ptr<T>(allocFunc(),[deleteFunc](T *p){ deleteFunc(p); });

}

};
}