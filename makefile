# Projects Options
CFLAGS = -Wall -Wextra -Weffc++ -Wunreachable-code -std=c++17
CFLAGS_DEBUG = -g
CFLAGS_NEEDED = -DIL_STD -lilocplex -lconcert -lcplex -lm -lpthread -ldl -m64
CFLAGS_RISKY = -O3 -funsafe-math-optimizations -funsafe-loop-optimizations -floop-nest-optimize -ffinite-math-only #Some few risky optimisation flags you can put in addons.
CFLAGS_ADDON = 
EXEC_NAME = exe

SYSTEM = x86-64_linux
LIBFORMAT = static_pic
CPLEXDIR = /home/local.isima.fr/tquvo/CPLEX_Studio1210/cplex
CONCERTDIR = /home/local.isima.fr/tquvo/CPLEX_Studio1210/concert
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
LEMONLIBDIR = /home/local.isima.fr/tquvo/lemon
INCLUDES = -I headers -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include -I$(LEMONLIBDIR)/include
LIBS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -m64 -lpthread -ldl -lstdc++fs -L$(LEMONLIBDIR)/lib -lemon

# CONFIGURATION
CC = g++
EXT = cpp
SRC_PATH1=src
SRC_PATH2=
OBJ_PATH=Objs

# Objects
SRC_FILES_REQUIRED = $(wildcard $(SRC_PATH1)/*.$(EXT) $(SRC_PATH2)/*.$(EXT))
OBJ_FILES_ALL = $(addprefix $(OBJ_PATH)/,$(notdir $(SRC_FILES_REQUIRED:.$(EXT)=.o)))
OBJ_FILES_ALL_DEBUG = $(OBJ_FILES_ALL:.o=.og)

OBJ_FILES_REQUIRED = $(filter-out $(OBJ_PATH)/main.o,$(filter-out $(OBJ_PATH)/test.o,$(filter-out $(OBJ_PATH)/catch.o,$(OBJ_FILES_ALL))))
OBJ_FILES_REQUIRED_DEBUG = $(OBJ_FILES_REQUIRED:.o=.og)

OBJ_FILES = $(OBJ_FILES_REQUIRED) $(OBJ_PATH)/main.o
OBJ_FILES_DEBUG = $(OBJ_FILES:.o=.og)
OBJ_FILES_TEST = $(OBJ_FILES_REQUIRED) $(OBJ_PATH)/test.o $(OBJ_PATH)/catch.o
OBJ_FILES_TEST_DEBUG = $(OBJ_FILES_TEST:.o=.og)
OBJ_FILES_MOSTLY = $(EXEC_NAME) $(OBJ_PATH)/VRP.o $(OBJ_PATH)/VRP.og $(OBJ_PATH)/test.o $(OBJ_PATH)/main.o $(OBJ_PATH)/test.og $(OBJ_PATH)/main.og

# Colors
NO_COLOR = \033[0m
OK_COLOR = \033[32;01m
ERROR_COLOR = \033[31;01m
WARN_COLOR = \033[33;01m

OK_STRING = $(OK_COLOR)[OK]$(NO_COLOR)
ERROR_STRING = $(ERROR_COLOR)[ERROR]$(NO_COLOR)
WARN_STRING = $(WARN_COLOR)[WARNING]$(NO_COLOR)

# Compilation Options
all : $(EXEC_NAME)
	@rm -f $(OBJ_FILES_DEBUG) test.og

help :
	@echo "Actual configuration:\n\tCompiler: $(CC) for $(EXT) files.\n\tSearching source files in $(SRC_PATH1)/ and $(SRC_PATH2)/ directories.\n\tObject files will be created in $(OBJ_PATH)/ directory.\n"
	@echo "make            :\tmain compilation (main.*)"
	@echo "make reduce     :\tmain compilation (main.*) and executable size reduction -> not debugable"
	@echo "make debug      :\tmain compilation (main.*) with debug flag"
	@echo "make test       :\tunit test (test.*) compilation"
	@echo "make testdebug  :\tunit test (test.*) compilation"
	@echo "make links      :\tlink step only"
	@echo "make clean      :\tclean all files"
	@echo "make mostlyclean:\tclean selected files"

reduce : $(EXEC_NAME)
	@rm -f $(OBJ_FILES_DEBUG) test.og
	@strip $(EXEC_NAME) -s -R .comment -R .gnu.version --strip-unneeded
	@echo "$(OK_STRING) $(EXEC_NAME) size reduced !"
	@echo "$(WARN_STRING) You will not be able to debug $(EXEC_NAME) anymore !"

debug : $(OBJ_FILES_DEBUG)
	@rm -f $(EXEC_NAME) $(OBJ_FILES) test.o
	@echo -n "compiling $(EXEC_NAME)...\n\033[1F"
	@$(CC) $(INCLUDES) -o $(EXEC_NAME) $(OBJ_FILES_DEBUG) $(LIBS)
	@echo "\033[2K$(OK_STRING) $(EXEC_NAME) compiled !"

test : $(OBJ_FILES_TEST)
	@rm -f $(EXEC_NAME) $(OBJ_FILES_DEBUG) test.og
	@echo -n "compiling $(EXEC_NAME)...\n\033[1F"
	@$(CC) $(INCLUDES) -o $(EXEC_NAME) $(OBJ_FILES_TEST) $(LIBS)
	@echo "\033[2K$(OK_STRING) $(EXEC_NAME) compiled !"

testdebug : $(OBJ_FILES_TEST_DEBUG)
	@rm -f $(EXEC_NAME) $(OBJ_FILES) test.o
	@echo -n "compiling $(EXEC_NAME)...\n\033[1F"
	@$(CC) $(INCLUDES) -o $(EXEC_NAME) $(OBJ_FILES_TEST_DEBUG) $(LIBS)
	@echo "\033[2K$(OK_STRING) $(EXEC_NAME) compiled !"

links :
	@echo -n "compiling $(EXEC_NAME)...\n\033[1F"
	@$(CC) $(INCLUDES) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)
	@echo "\033[2K$(OK_STRING) $(EXEC_NAME) compiled !"

clean :
	@rm -f $(OBJ_FILES_ALL) $(OBJ_FILES_ALL_DEBUG)
	@echo "$(OK_STRING) all files cleaned !"

mostlyclean :
	@rm -f $(OBJ_FILES_MOSTLY)
	@echo "$(OK_STRING) $(OBJ_FILES_MOSTLY) cleaned !"

$(EXEC_NAME) : $(OBJ_FILES)
	@echo -n "compiling $@...\n\033[1F"
	@$(CC) $(INCLUDES) -o $(EXEC_NAME) $(OBJ_FILES) $(LIBS)
	@echo "\033[2K$(OK_STRING) $@ compiled !"

$(OBJ_PATH)/%.o : $(SRC_PATH1)/%.$(EXT) | $(OBJ_PATH)
	@echo -n "compiling $@...\n\033[1F"
	@$(CC) $(INCLUDES) -o $@ -c $^ $(CFLAGS) $(CFLAGS_NEEDED) $(CFLAGS_ADDON)
	@echo "\033[2K$(OK_STRING) $@ compiled !"

$(OBJ_PATH)/%.og : $(SRC_PATH2)/%.$(EXT) | $(OBJ_PATH)
	@echo -n "compiling $@...\n\033[1F"
	@$(CC) $(INCLUDES) -o $@ -c $^ $(CFLAGS) $(CFLAGS_NEEDED) $(CFLAGS_DEBUG)
	@echo "\033[2K$(OK_STRING) $@ compiled !"

$(OBJ_PATH)/%.o : $(SRC_PATH2)/%.$(EXT) | $(OBJ_PATH)
	@echo -n "compiling $@...\n\033[1F"
	@$(CC) $(INCLUDES) -o $@ -c $^ $(CFLAGS) $(CFLAGS_NEEDED) $(CFLAGS_ADDON)
	@echo "\033[2K$(OK_STRING) $@ compiled !"

$(OBJ_PATH)/%.og : $(SRC_PATH1)/%.$(EXT) | $(OBJ_PATH)
	@echo -n "compiling $@...\n\033[1F"
	@$(CC) $(INCLUDES) -o $@ -c $^ $(CFLAGS) $(CFLAGS_NEEDED) $(CFLAGS_DEBUG)
	@echo "\033[2K$(OK_STRING) $@ compiled !"

$(OBJ_PATH) :
	mkdir $(OBJ_PATH)
