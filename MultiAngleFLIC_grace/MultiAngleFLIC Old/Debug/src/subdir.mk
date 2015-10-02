################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/MultiAngleFLIC.cpp 

OBJS += \
./src/MultiAngleFLIC.o 

CPP_DEPS += \
./src/MultiAngleFLIC.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel IA-32 C++ Compiler '
	icpc -g -I/opt/intel/Compiler/11.1/072/mkl/include -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


