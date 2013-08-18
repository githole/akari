#ifndef _SETTING_H_
#define _SETTING_H_

#include <iostream>
#include <string>
#include <fstream>
#include <map>


class Setting {
private:
	std::map<std::string, std::string> value_table_;
public:
	Setting(const std::string &filename) {
		std::ifstream ifs(filename);

		if (!ifs)
			return;

		while (!ifs.eof()) {
			std::string key;
			std::string value;
			ifs >> key >> value;

			value_table_.insert(std::make_pair(key, value));
		}
	}

	float float_value(const std::string &key, const float default_value = 0.0f) {
		if (value_table_.find(key) == value_table_.end())
			return default_value;
		return (float)atof(value_table_[key].c_str());
	}
	
	int int_value(const std::string &key, const int default_value = 0) {
		if (value_table_.find(key) == value_table_.end())
			return default_value;
		return atoi(value_table_[key].c_str());
	}
	
	std::string string_value(const std::string &key, const std::string &default_value = "") {
		if (value_table_.find(key) == value_table_.end())
			return default_value;
		return value_table_[key];
	}
};

#endif