/** Logger
 *
 * @author Markus Mirz <mmirz@eonerc.rwth-aachen.de>
 * @copyright 2017, Institute for Automation of Complex Power Systems, EONERC
 * @license GNU General Public License (version 3)
 *
 * DPsim
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *********************************************************************************/

#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;


#include "Logger.h"

using namespace DPsim;

std::ostringstream Logger::nullStream;

Logger::Logger() : mLogFile() {
	mLogLevel = LogLevel::NONE;
	mLogFile.setstate(std::ios_base::badbit);
}

Logger::Logger(String filename, LogLevel level) : mLogLevel(level) {
	fs::path p = filename;

	if (p.has_parent_path() && !fs::exists(p.parent_path()))
		fs::create_directory(p.parent_path());

	mLogFile = std::ofstream(filename);
	if (!mLogFile.is_open()) {
		std::cerr << "Cannot open log file " << filename << std::endl;
		mLogLevel = LogLevel::NONE;
	}
}

Logger::~Logger() {
	if (mLogFile.is_open())
		mLogFile.close();
}

std::ostream& Logger::Log(LogLevel level) {
	if (level > mLogLevel) {
		return getNullStream();
	}

	switch (level) {
		case LogLevel::INFO:
			mLogFile << "INFO: ";
			break;
		case LogLevel::WARN:
			mLogFile << "WARN: ";
			break;
		case LogLevel::ERROR:
			mLogFile << "ERROR: ";
			break;
		case LogLevel::NONE:
			return getNullStream();
			break;
	}
	return mLogFile;
}

void Logger::Log(LogLevel level, String str) {
	if (level > mLogLevel) {
		return;
	}

	switch (level) {
	case LogLevel::INFO:
		mLogFile << "INFO: " << str << std::endl;
		break;
	case LogLevel::WARN:
		mLogFile << "WARN: " << str << std::endl;
		break;
	case LogLevel::ERROR:
		mLogFile << "ERROR: " << str << std::endl;
		break;
	case LogLevel::NONE:
		return;
		break;
	}
}

void Logger::LogMatrix(LogLevel level, Matrix& data) {
	if (level > mLogLevel) {
		return;
	}

	mLogFile << data << std::endl;
}

void Logger::LogMatrix(LogLevel level, const Matrix& data) {
	if (level > mLogLevel) {
		return;
	}

	mLogFile << data << std::endl;
}

void Logger::LogDataLine(Real time, Matrix& data) {
	mLogFile << std::scientific << time;

	for (Int i = 0; i < data.rows(); i++)
		mLogFile << ", " << data(i, 0);

	mLogFile << std::endl;
}

void Logger::LogDataLine(Real time, Real data) {
	mLogFile << std::scientific << time;
	mLogFile << ", " << data;
	mLogFile << std::endl;
}
std::ostream& Logger::getNullStream() {
	if (nullStream.good())
		nullStream.setstate(std::ios_base::badbit);
	return nullStream;
}
