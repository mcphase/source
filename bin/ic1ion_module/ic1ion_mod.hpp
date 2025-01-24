#pragma once
#include <iostream>
#include "singleion_module.hpp"
class ic1ion_mod : public singleion_module
{
public:
ic1ion_mod() = default;
~ic1ion_mod() = default;
void Icalc() override;
};
