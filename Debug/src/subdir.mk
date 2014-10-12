################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/fdct3d_forward.o \
../src/fdct3d_inverse.o \
../src/fdct3d_param.o \
../src/test.o 

CPP_SRCS += \
../src/fdct3d_forward.cpp \
../src/fdct3d_inverse.cpp \
../src/fdct3d_param.cpp \
../src/test.cpp 

OBJS += \
./src/fdct3d_forward.o \
./src/fdct3d_inverse.o \
./src/fdct3d_param.o \
./src/test.o 

CPP_DEPS += \
./src/fdct3d_forward.d \
./src/fdct3d_inverse.d \
./src/fdct3d_param.d \
./src/test.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -Ifftw3 -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


