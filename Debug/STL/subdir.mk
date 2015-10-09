################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../STL/Always.cpp \
../STL/Atom.cpp \
../STL/Conjunction.cpp \
../STL/Disjunction.cpp \
../STL/Eventually.cpp \
../STL/STL.cpp \
../STL/Until.cpp 

OBJS += \
./STL/Always.o \
./STL/Atom.o \
./STL/Conjunction.o \
./STL/Disjunction.o \
./STL/Eventually.o \
./STL/STL.o \
./STL/Until.o 

CPP_DEPS += \
./STL/Always.d \
./STL/Atom.d \
./STL/Conjunction.d \
./STL/Disjunction.d \
./STL/Eventually.d \
./STL/STL.d \
./STL/Until.d 


# Each subdirectory must supply rules for building sources it contributes
STL/%.o: ../STL/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/usr/local/include/CGAL -I/home/dreossi/eigen-eigen-bdd17ee3b1b3 -O0 -g3 -Wall -c -fmessage-length=0 --std=c++0x -frounding-math -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


