################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../BaseConverter.cpp \
../Box.cpp \
../Bundle.cpp \
../CoeffMap.cpp \
../DiscreteDynamicalSystem.cpp \
../DynamicalSystem.cpp \
../LinearSystem.cpp \
../LinearSystemSet.cpp \
../Model.cpp \
../MultiParallelotope.cpp \
../MultiReacher.cpp \
../Parallelotope.cpp \
../ParameterSynthesizer.cpp \
../Polyhedron.cpp \
../VarsGenerator.cpp \
../main_ebola.cpp \
../main_influenza_multi.cpp \
../main_old_sir.cpp \
../main_osc_multi.cpp \
../main_quad_multi.cpp \
../main_sars_multi.cpp \
../main_si.cpp \
../main_sir.cpp \
../main_sir_multi.cpp \
../main_thesis.cpp 

OBJS += \
./BaseConverter.o \
./Box.o \
./Bundle.o \
./CoeffMap.o \
./DiscreteDynamicalSystem.o \
./DynamicalSystem.o \
./LinearSystem.o \
./LinearSystemSet.o \
./Model.o \
./MultiParallelotope.o \
./MultiReacher.o \
./Parallelotope.o \
./ParameterSynthesizer.o \
./Polyhedron.o \
./VarsGenerator.o \
./main_ebola.o \
./main_influenza_multi.o \
./main_old_sir.o \
./main_osc_multi.o \
./main_quad_multi.o \
./main_sars_multi.o \
./main_si.o \
./main_sir.o \
./main_sir_multi.o \
./main_thesis.o 

CPP_DEPS += \
./BaseConverter.d \
./Box.d \
./Bundle.d \
./CoeffMap.d \
./DiscreteDynamicalSystem.d \
./DynamicalSystem.d \
./LinearSystem.d \
./LinearSystemSet.d \
./Model.d \
./MultiParallelotope.d \
./MultiReacher.d \
./Parallelotope.d \
./ParameterSynthesizer.d \
./Polyhedron.d \
./VarsGenerator.d \
./main_ebola.d \
./main_influenza_multi.d \
./main_old_sir.d \
./main_osc_multi.d \
./main_quad_multi.d \
./main_sars_multi.d \
./main_si.d \
./main_sir.d \
./main_sir_multi.d \
./main_thesis.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/usr/local/include/CGAL -I/home/dreossi/eigen-eigen-bdd17ee3b1b3 -O0 -g3 -Wall -c -fmessage-length=0 --std=c++0x -frounding-math -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


