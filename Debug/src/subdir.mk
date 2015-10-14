################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Always.cpp \
../src/Atom.cpp \
../src/BaseConverter.cpp \
../src/Bundle.cpp \
../src/Conjunction.cpp \
../src/Disjunction.cpp \
../src/Eventually.cpp \
../src/LinearSystem.cpp \
../src/LinearSystemSet.cpp \
../src/Parallelotope.cpp \
../src/STL.cpp \
../src/Sapo.cpp \
../src/Until.cpp \
../src/VarsGenerator.cpp \
../src/main_osc_multi.cpp \
../src/main_quad_multi2.cpp \
../src/main_sir_multi.cpp 

OBJS += \
./src/Always.o \
./src/Atom.o \
./src/BaseConverter.o \
./src/Bundle.o \
./src/Conjunction.o \
./src/Disjunction.o \
./src/Eventually.o \
./src/LinearSystem.o \
./src/LinearSystemSet.o \
./src/Parallelotope.o \
./src/STL.o \
./src/Sapo.o \
./src/Until.o \
./src/VarsGenerator.o \
./src/main_osc_multi.o \
./src/main_quad_multi2.o \
./src/main_sir_multi.o 

CPP_DEPS += \
./src/Always.d \
./src/Atom.d \
./src/BaseConverter.d \
./src/Bundle.d \
./src/Conjunction.d \
./src/Disjunction.d \
./src/Eventually.d \
./src/LinearSystem.d \
./src/LinearSystemSet.d \
./src/Parallelotope.d \
./src/STL.d \
./src/Sapo.d \
./src/Until.d \
./src/VarsGenerator.d \
./src/main_osc_multi.d \
./src/main_quad_multi2.d \
./src/main_sir_multi.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/usr/local/include/CGAL -I../include/ -I/home/dreossi/eigen-eigen-bdd17ee3b1b3 -O0 -g3 -Wall -c -fmessage-length=0 --std=c++0x -frounding-math -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


