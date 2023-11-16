#include <BAT/BCModel.h>
#include <BAT/BCMath.h>
#include <BAT/BCMTF.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>

#include <TFile.h>
#include <TMath.h>

#include <string>
#include <vector>
#include "../settings/json.hpp"
using namespace nlohmann;
#include <sys/stat.h>

#pragma once

        bool fileExists(const std::string& filePath); // check if a file exists, given the path
        void write_to_json(std::string filePath, std::string name_fit, ordered_json foutput); // write output to json file
	std::vector<int> CheckPosterior(std::string fit_name, std::string E0, int bkg); // check posterior pdfs
