/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2015, The University of Texas and Bryan Marker

    DxTer is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DxTer is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DxTer.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <fstream>

#include "logging.h"
#include "driverUtils.h"

using namespace std;

ofstream* logStream;

string log_dir = "logs";

string LogLevelToString(LogLevel level) {
  switch(level) {
  case(ALL):
    return "ALL";
  case(INFO):
    return "INFO";
  case(DEBUG):
    return "DEBUG";
  case(FATAL):
    return "FATAL";
  default:
    LOG_FAIL("Bad case in LogLevelToString");
  }
}

void OpenLog(string logName) {
  string logFileName = log_dir + "/" + logName + "-" + DateAndTimeString() + ".dxt_log";
  logStream = new ofstream(logFileName, ofstream::out);
  LOG_A("Started log in file " + logFileName);
}

void CloseLog() {
  LOG_A("End of log");
  logStream->close();
  delete logStream;
}

void LogMessage(LogLevel level, string message) {
  string msg = "";
  msg += LogLevelToString(level) + "\t";
  msg += DateAndTimeString();
  msg += ": " + message;
  *logStream << msg << endl;
}

void LogFailMessageAndThrow(string message) {
  LOG_F(message);
  END_LOG();
  throw;
}
