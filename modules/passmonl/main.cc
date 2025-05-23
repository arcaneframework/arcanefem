#include <arcane/launcher/ArcaneLauncher.h>

using namespace Arcane;

int
main(int argc,char* argv[])
{
    ArcaneLauncher::init(CommandLineArguments(&argc,&argv));
    auto& app_build_info = ArcaneLauncher::applicationBuildInfo();
    app_build_info.setCodeName("Passmonl");
    app_build_info.setCodeVersion(VersionInfo(1,0,0));
    return ArcaneLauncher::run();
}
