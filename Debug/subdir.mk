################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../BaseConverter.cpp \
../Box.cpp \
../DiscreteDynamicalSystem.cpp \
../DynamicalSystem.cpp \
../LinearSystem.cpp \
../LinearSystemSet.cpp \
../Model.cpp \
../Parallelotope.cpp \
../ParameterSynthesizer.cpp \
../Polyhedron.cpp \
../main_ebola.cpp \
../main_si.cpp 

OBJS += \
./BaseConverter.o \
./Box.o \
./DiscreteDynamicalSystem.o \
./DynamicalSystem.o \
./LinearSystem.o \
./LinearSystemSet.o \
./Model.o \
./Parallelotope.o \
./ParameterSynthesizer.o \
./Polyhedron.o \
./main_ebola.o \
./main_si.o 

CPP_DEPS += \
./BaseConverter.d \
./Box.d \
./DiscreteDynamicalSystem.d \
./DynamicalSystem.d \
./LinearSystem.d \
./LinearSystemSet.d \
./Model.d \
./Parallelotope.d \
./ParameterSynthesizer.d \
./Polyhedron.d \
./main_ebola.d \
./main_si.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/usr/local/include/CGAL -O0 -g3 -Wall -c -fmessage-length=0 --std=c++0x -frounding-math -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


