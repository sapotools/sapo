################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Atom.cpp \
../BaseConverter.cpp \
../Box.cpp \
../Conjunction.cpp \
../DiscreteDynamicalSystem.cpp \
../Disjunction.cpp \
../DynamicalSystem.cpp \
../Ebola.cpp \
../LinearSystem.cpp \
../LinearSystemSet.cpp \
../Model.cpp \
../Parallelotope.cpp \
../ParameterSynthesizer.cpp \
../Polyhedron.cpp \
../STL.cpp \
../Until.cpp \
../main.cpp 

OBJS += \
./Atom.o \
./BaseConverter.o \
./Box.o \
./Conjunction.o \
./DiscreteDynamicalSystem.o \
./Disjunction.o \
./DynamicalSystem.o \
./Ebola.o \
./LinearSystem.o \
./LinearSystemSet.o \
./Model.o \
./Parallelotope.o \
./ParameterSynthesizer.o \
./Polyhedron.o \
./STL.o \
./Until.o \
./main.o 

CPP_DEPS += \
./Atom.d \
./BaseConverter.d \
./Box.d \
./Conjunction.d \
./DiscreteDynamicalSystem.d \
./Disjunction.d \
./DynamicalSystem.d \
./Ebola.d \
./LinearSystem.d \
./LinearSystemSet.d \
./Model.d \
./Parallelotope.d \
./ParameterSynthesizer.d \
./Polyhedron.d \
./STL.d \
./Until.d \
./main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


