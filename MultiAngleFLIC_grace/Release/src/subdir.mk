################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/AuxiliaryFXNs.cpp \
../src/MultiAngleFLIC.cpp \
../src/Surface.cpp \
../src/spectra.cpp 

OBJS += \
./src/AuxiliaryFXNs.o \
./src/MultiAngleFLIC.o \
./src/Surface.o \
./src/spectra.o 

CPP_DEPS += \
./src/AuxiliaryFXNs.d \
./src/MultiAngleFLIC.d \
./src/Surface.d \
./src/spectra.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	icpc -inline-level=2 -I/opt/intel/Compiler/11.1/073/mkl/include -I/usr/include/ImageMagick -xSSE4.1 -O2 -std=c99 -openmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


