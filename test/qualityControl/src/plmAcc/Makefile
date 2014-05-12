CC 	= gcc
CLFAGS	:= -Wall -fPIC
LINK	:= -shared -Wl
INC	:= -I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux
TARGET	= jniWrapper 

OBJS	= 

ifeq ($(dbg),1)
	CFLAGS += -g
else
	CFLAGS += -O3
endif


all:
	$(CC) $(CFLAGS) -o lib.so $(INC) %.c

clean:
	rm *.o *.so
