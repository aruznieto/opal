//
// Created by anrunie on 8/11/23.
//

#ifndef OPALDEV_CONFIGURATION_H
#define OPALDEV_CONFIGURATION_H

#include "json.hpp"
#include <vector>
#include <string>
#include <utility>
using json=nlohmann::json;


//Read configuration from a JSON file.
//Format is flexible, the user can read any key  defined in the config file passed as parameter
namespace opal {
    class ConfigurationParser {
        public:
            ConfigurationParser(const std::string& configFile);
            nlohmann::json getKey(const std::string& key) const;
            nlohmann::json getKey(const nlohmann::json& node, const std::string& key) const;
        private:
            nlohmann::json json_;
    };
}

#endif //OPALDEV_CONFIGURATION_H
