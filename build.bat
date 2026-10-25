@echo off
rem Build with MinGW g++ (adjust the compiler name if needed)
g++ -O2 -std=c++17 -Wall -Wextra -o alpha_wrap.exe src\alpha_wrap.cpp
if errorlevel 1 (
  echo Build failed.
  exit /b 1
)
echo Built alpha_wrap.exe
