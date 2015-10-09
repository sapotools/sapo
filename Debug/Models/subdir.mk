################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Models/Converge.cpp \
../Models/Ebola.cpp \
../Models/SI.cpp \
../Models/SIR.cpp 

OBJS += \
./Models/Converge.o \
./Models/Ebola.o \
./Models/SI.o \
./Models/SIR.o 

CPP_DEPS += \
./Models/Converge.d \
./Models/Ebola.d \
./Models/SI.d \
./Models/SIR.d 


# Each subdirectory must supply rules for building sources it contributes
Models/%.o: ../Models/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/usr/local/include/CGAL -I/home/dreossi/eigen-eigen-bdd17ee3b1b3 -O0 -g3 -Wall -c -fmessage-length=0 --std=c++0x -frounding-math -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


