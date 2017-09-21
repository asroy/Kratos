#pragma once
#include<cxxabi.h>
#include<typeinfo>

namespace DataUtility
{

template<typename T>
void type_name(T u)
{
  int     status;
  char   *realname;

  const std::type_info  &ti = typeid(u);

  realname = abi::__cxa_demangle(ti.name(), 0, 0, &status);
  std::cout << ti.name() << "\t=> " << realname << "\t: " << status << '\n';
  free(realname);
}

}