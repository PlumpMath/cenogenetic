#adapted from Learn C the Hard Way

ifeq "$(CXX)" "g++"
	OPTFLAGS+=-Wno-parentheses
endif
ifeq "$(CXX)" "clang++"
	OPTFLAGS+=-Wno-dangling-else
endif
CFLAGS=-g -O2 -Wall -Wextra -Isrc -rdynamic -DNDEBUG $(OPTFLAGS)
CXXFLAGS= -std=c++11 $(CFLAGS)
LIBS=-ldl $(OPTLIBS) PREFIX?=/usr/local

SOURCES=$(wildcard src/**/*.cc src/*.cc)
OBJECTS=$(patsubst %.cc,%.o,$(SOURCES))

TEST_SRC=$(wildcard tests/*_tests.cc)
TESTS=$(patsubst %.cc,%,$(TEST_SRC))

TARGET=build/libcenogenetic.a
SO_TARGET=$(patsubst %.a,%.so,$(TARGET))

all: pretty tags libs tests

libs: $(TARGET) $(SO_TARGET)

dev: CFLAGS=-g -Wall -Isrc -Wall -Wextra $(OPTFLAGS)
dev: CXXFLAGS= -std=c++11 $(CFLAGS)
dev: libs tests

$(TARGET): CFLAGS += -fPIC
$(TARGET): build $(OBJECTS)
	ar rcs $@ $(OBJECTS)
	ranlib $@

$(SO_TARGET): $(TARGET) $(OBJECTS)
	$(CC) -shared -o $@ $(OBJECTS)

build:
	@mkdir -p build
	@mkdir -p bin

.PHONY: tests
#tests: CFLAGS += $(TARGET)
tests: CFLAGS += $(OBJECTS)
tests: $(TESTS)
	sh ./tests/runtests.sh

valgrind:
	VALGRIND="valgrind --log-file=/tmp/valgrind-%p.log" $(MAKE)

clean:
	rm -rf build $(OBJECTS) $(TESTS)
	rm -f tests/tests.log
	find . -name "*.gc*" -exec rm {} \;
	rm -rf `find . -name "*.dSYM" -print`
	find . -name "*~" -exec rm {} \;
	find . -name "*.orig" -exec rm {} \;
	rm -f tags

install: all
	install -d $(DESTDIR)/$(PREFIX)/lib/
	install $(TARGET) $(DESTDIR)/$(PREFIX)/lib/

.PHONY: pretty
pretty:
	if which astyle >/dev/null ; then astyle --options=./astylerc $(TEST_SRC) $(SOURCES) src/*.h ; fi
	#if which vim >/dev/null ; then vim -c ':bufdo normal gg=G' -c ':xall' $(TEST_SRC) $(SOURCES) src/*.h ; fi

.PHONY: tags
tags:
	if which ctags >/dev/null ; then ctags -R; fi

BADFUNCS='[^_.>a-zA-Z0-9](str(n?cpy|n?cat|xfrm|n?dup|str|pbrk|tok|_)|stpn?cpy|a?sn?printf|byte_)'
check:
	@echo Files with potentially dangerous functions.
	@egrep $(BADFUNCS) $(SOURCES) || true
