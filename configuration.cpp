//
// Created by anrunie on 8/11/23.
//

#include "configuration.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>

namespace opal {

    ConfigurationParser::ConfigurationParser(const std::string& configFile) {
        std::ifstream file(configFile);
        if (file.is_open()) {
            try {
                std::cout << "Reading JSON: " << configFile << '\n';
                file >> json_;
            } catch (const nlohmann::json::parse_error& e) {
                std::cerr << "JSON parse error: " << e.what() << '\n';
            }
        } else {
            std::cerr << "Unable to open config file: " << configFile << '\n';
        }
    }

    // Get a value by key from the JSON file at the root level or from a specified node
    nlohmann::json ConfigurationParser::getKey(const std::string& key) const {
        if (json_.contains(key)) {
            return getKey(json_, key);
        } else {
            std::cerr << "Key not found: " << key << '\n';
            return nullptr; // Or an appropriate default value or error indication
        }
    }

    nlohmann::json ConfigurationParser::getKey(const nlohmann::json& node, const std::string& key) const {
        if (node.contains(key)) {
            return node.at(key);
        } else {
            std::cerr << "Key not found: " << key << '\n';
            return nullptr; // Or an appropriate default value or error indication
        }
    }

}