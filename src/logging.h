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

#include <string>

#ifndef LOGGING_H
#define LOGGING_H

using namespace std;

enum LogLevel { ALL, INFO, DEBUGVAL, FATAL };

#define LOG_START(name) OpenLog((name))
#define LOG_END() CloseLog()
#define LOG_A(msg) LogMessage(ALL, (msg))
#define LOG_I(msg) LogMessage(INFO, (msg))
#define LOG_D(msg) LogMessage(DEBUGVAL, (msg))
#define LOG_F(msg) LogMessage(FATAL, (msg))
#define LOG_FAIL(failMsg) LogFailMessageAndThrow((failMsg));

void OpenLog(string logName);
void CloseLog();
void LogMessage(LogLevel level, string message);
void LogFailMessageAndThrow(string message);

#endif // LOGGING_H
