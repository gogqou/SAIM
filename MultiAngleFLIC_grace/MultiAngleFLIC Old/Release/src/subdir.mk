################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/AuxiliarFXNs.cpp \
../src/MultiAngleFLIC.cpp 

OBJS += \
./src/AuxiliarFXNs.o \
./src/MultiAngleFLIC.o 

CPP_DEPS += \
./src/AuxiliarFXNs.d \
./src/MultiAngleFLIC.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel IA-32 C++ Compiler '
	icpc -inline-level=2 -mtune=core2 -I/opt/intel/Compiler/11.1/072/mkl/include -I/home/matt/workspace/MultiAngleFLIC/src -I/usr/include/ImageMagick -O2 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


