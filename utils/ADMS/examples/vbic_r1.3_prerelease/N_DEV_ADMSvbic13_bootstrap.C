#include <Xyce_config.h>
#include <N_DEV_ADMSvbic13.h>
#include <N_DEV_ADMSvbic13_4t.h>

struct Bootstrap
{
  Bootstrap()
  {
    Xyce::Device::ADMSvbic13::registerDevice();
    Xyce::Device::ADMSvbic13_4t::registerDevice();
  }
};

Bootstrap s_bootstrap;
