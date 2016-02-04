#include "Config/ConfigurationLoader.h"
#include <string.h>

#include "Config/TemplateWanderConfigurationLoader.h"
#include "Config/TemplateBoidsConfigurationLoader.h"
#include "Config/TemplateMedeaConfigurationLoader.h"
#include "Config/TemplateRandomwalkConfigurationLoader.h"
#include "Config/IncrementConfigurationLoader.h"
#include "Config/TestNNConfigurationLoader.h"
//###DO-NOT-DELETE-THIS-LINE###TAG:INCLUDE###//


ConfigurationLoader::ConfigurationLoader()
{
	//nothing to do
}

ConfigurationLoader::~ConfigurationLoader()
{
	//nothing to do
}

ConfigurationLoader* ConfigurationLoader::make_ConfigurationLoader (std::string configurationLoaderObjectName)
{
	if (0)
	{
		// >>> Never reached
	}
#if defined PRJ_TEMPLATEWANDER || !defined MODULAR
	else if (configurationLoaderObjectName == "TemplateWanderConfigurationLoader" )
	{
		return new TemplateWanderConfigurationLoader();
	}
#endif
#if defined PRJ_TEMPLATEBOIDS || !defined MODULAR
	else if (configurationLoaderObjectName == "TemplateBoidsConfigurationLoader" )
	{
		return new TemplateBoidsConfigurationLoader();
	}
#endif
#if defined PRJ_TEMPLATEMEDEA || !defined MODULAR
	else if (configurationLoaderObjectName == "TemplateMedeaConfigurationLoader" )
	{
		return new TemplateMedeaConfigurationLoader();
	}
#endif
#if defined PRJ_TEMPLATERANDOMWALK || !defined MODULAR
	else if (configurationLoaderObjectName == "TemplateRandomwalkConfigurationLoader" )
	{
		return new TemplateRandomwalkConfigurationLoader();
	}
#endif
#if defined PRJ_INCREMENT || !defined MODULAR
	else if (configurationLoaderObjectName == "IncrementConfigurationLoader" )
	{
		return new IncrementConfigurationLoader();
	}
#endif
#if defined PRJ_TESTNN || !defined MODULAR
	else if (configurationLoaderObjectName == "TestNNConfigurationLoader" )
	{
		return new TestNNConfigurationLoader();
	}
#endif
    //###DO-NOT-DELETE-THIS-LINE###TAG:SWITCH###//
	else
	{
		return NULL;
	}

}
