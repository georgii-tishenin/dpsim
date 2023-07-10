/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <memory>
#include <iomanip>

#include <spdlog/sinks/null_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <dpsim-models/Filesystem.h>
#include <dpsim-models/Logger.h>

using namespace CPS;

spdlog::sink_ptr Logger::mCliSink;
spdlog::sink_ptr Logger::mComponentFileSink;

void Logger::setLogLevel(std::shared_ptr<spdlog::logger> logger, Logger::Level level) {
	logger->set_level(level);
}

void Logger::setLogPattern(std::shared_ptr<spdlog::logger> logger, std::string pattern) {
	logger->set_pattern(pattern);
}

// #### to string methods ####
String Logger::matrixToString(const Matrix& mat) {
	std::stringstream ss;
	ss << std::scientific << "\n" << mat;
	return ss.str();
}

String Logger::matrixCompToString(const MatrixComp& mat) {
	std::stringstream ss;
	ss << std::scientific << "\n" << mat;
	return ss.str();
}

String Logger::sparseMatrixToString(const SparseMatrix& mat) {
	return matrixToString(Matrix(mat));
}

String Logger::sparseMatrixCompToString(const SparseMatrixComp& mat) {
	return matrixCompToString(MatrixComp(mat));
}

String Logger::phasorMatrixToString(const MatrixComp& mat) {
	std::stringstream ss;
	ss << std::scientific << Math::abs(mat) << "\n\n" << Math::phase(mat);
	return ss.str();
}

String Logger::phasorToString(const Complex& num) {
	std::stringstream ss;
	ss << std::defaultfloat << Math::abs(num) << "<" << Math::phaseDeg(num);
	return ss.str();
}

String Logger::complexToString(const Complex& num) {
	std::stringstream ss;
	ss << std::defaultfloat << num.real() << "+j" << num.imag();
	return ss.str();
}

String Logger::realToString(const Real& num) {
	std::stringstream ss;
	ss << std::defaultfloat << num;
	return ss.str();
}

String Logger::prefix() {
	char const *p = getenv("CPS_LOG_PREFIX");

	return p ? p : "";
}

String Logger::logDir() {
	char const *p = getenv("CPS_LOG_DIR");

	return p ? p : "logs";
}

/// Set env variable CPS_LOG_DIR and overwrite
void Logger::setLogDir(String path) {
#ifdef __linux__
	setenv("CPS_LOG_DIR", path.c_str(), 1);
#elif defined(_WIN32)
	String var = "CPS_LOG_DIR=" + path;
	_putenv(var.c_str());
#endif
}

String Logger::getCSVColumnNames(std::vector<String> names) {
	std::stringstream ss;
    ss << std::right << std::setw(14) << "time";
	for (auto const& name : names) {
		ss << ", " << std::right << std::setw(13) << name;
	}
	ss << '\n';

	return ss.str();
}

String Logger::getCSVLineFromData(Real time, Real data) {
	std::stringstream ss;
    ss << std::scientific << std::right << std::setw(14) << time;
	ss << ", " << std::right << std::setw(13) << data;
	ss << '\n';

	return ss.str();
}

String Logger::getCSVLineFromData(Real time, const Matrix& data) {
	std::stringstream ss;
    ss << std::scientific << std::right << std::setw(14) << time;
	for (Int i = 0; i < data.rows(); i++) {
		ss << ", " << std::right << std::setw(13) << data(i, 0);
	}
	ss << '\n';

	return ss.str();
}

String Logger::getCSVLineFromData(Real time, const MatrixComp& data) {
	std::stringstream ss;
    ss << std::scientific << std::right << std::setw(14) << time;
	for (Int i = 0; i < data.rows(); i++) {
		ss << ", " << std::right << std::setw(13) << data(i, 0);
	}
	ss << '\n';

	return ss.str();
}

Logger::Log Logger::get(LoggerType type, const std::string &name, Level filelevel, Level clilevel) {
	switch (type) {
		case LoggerType::SIMULATION: {
			if (Logger::Log logger = spdlog::get("simulation")) {
				return logger;
			} else return create(LoggerType::SIMULATION, "simulation", name, filelevel, clilevel);
		}
		case LoggerType::COMPONENT: {
			if (Logger::Log logger = spdlog::get(name)) {
				return logger;
			} else return create(LoggerType::COMPONENT, name, "components", filelevel, clilevel);
		}
		case LoggerType::DEBUG: {
			if (Logger::Log logger = spdlog::get(name)) {
				return logger;
			} else return create(LoggerType::DEBUG, name, name, filelevel, clilevel);
		}
		default: {
			// UNREACHABLE!
			return nullptr;
		}
	}
}

Logger::Log Logger::create(Logger::LoggerType type, const std::string &name, const std::string &fileName, Level filelevel, Level clilevel) {
	String logDir = Logger::logDir();
	String filepath = logDir + "/" + fileName + ".log";
	std::vector<spdlog::sink_ptr> sinks;
	Logger::Log ret;

	if (clilevel != Logger::Level::off) {
		// Why do we log into stderr instead of stdout?
		if (!Logger::mCliSink) {
			Logger::mCliSink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
			Logger::mCliSink->set_level(clilevel);
			Logger::mCliSink->set_pattern(fmt::format("{}[%T.%f %n %^%l%$] %v", CPS::Logger::prefix()));
		}
		sinks.push_back(Logger::mCliSink);
	}

	if (filelevel != Logger::Level::off) {
		spdlog::sink_ptr file_sink;

		if (type == LoggerType::COMPONENT) {
			// Use common file sink
			if (!Logger::mComponentFileSink) {
				Logger::mComponentFileSink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(filepath, true);
				Logger::mComponentFileSink->set_level(filelevel);
				// TODO: Use better prefix
				Logger::mComponentFileSink->set_pattern(prefix() + "[%l][%n] %v");
			}
			file_sink = Logger::mComponentFileSink;
		} else {
			// Create new file sink
			file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(filepath, true);
			file_sink->set_level(filelevel);
			// TODO: Use better prefix
			file_sink->set_pattern(prefix() + "[%L] %v");
		}
		sinks.push_back(file_sink);
	}

	if (filelevel == Logger::Level::off && clilevel == Logger::Level::off) {
		ret = spdlog::create<spdlog::sinks::null_sink_st>(name);
	} else {
		ret = std::make_shared<spdlog::logger>(name, begin(sinks), end(sinks));
		spdlog::register_logger(ret);
	}

	return ret;
}
