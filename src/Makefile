CC=gcc
FLAGS=-std=c11 -Wall -Werror -Wextra
LIBS_MACOS= -L /usr/local/lib -lcheck
LIBS_LINUX= -lcheck -lm -lpthread -lrt -lsubunit -D_GNU_SOURCE
 
OS = $(shell uname)
ifeq ($(OS), Darwin)
OPEN=open
LIBS=$(LIBS_MACOS)
# LEAK_CHECK = leaks --atExit -- ./test
endif
 
ifeq ($(OS), Linux)
OPEN=xdg-open
LIBS=$(LIBS_LINUX)
# LEAK_CHECK = valgrind --leak-check=full --show-leak-kinds=all -s ./test
endif

all: s21_math.a 

s21_math.o_with_gcov: s21_math.c
	$(CC) $(FLAGS) -c s21_math.c --coverage

s21_math.a_with_gcov: s21_math.o_with_gcov
	ar -rcs s21_math.a s21_math.o

test: clean s21_test.check s21_math.a_with_gcov
	checkmk clean_mode=1 s21_test.check > s21_test.c
	$(CC) s21_test.c s21_math.a -o s21_test $(LIBS) --coverage
	./s21_test
	rm -f s21_test.gcda s21_test.gcno s21_test.c

test_val: clean s21_test.check s21_math.a_with_gcov
	checkmk clean_mode=1 s21_test.check > s21_test.c
	$(CC) s21_test.c s21_math.a -o s21_test $(LIBS) --coverage
	valgrind ./s21_test > log.txt 2>&1
	grep -e "ERROR SUMMARY" log.txt | grep -v -e "ERROR SUMMARY: 0"
	rm -f s21_test.gcda s21_test.gcno s21_test.c

s21_math.a: s21_math.c
	$(CC) -c s21_math.c
	ar -rcs s21_math.a s21_math.o

gcov_report: test
	gcovr --html --html-details -o report.html
	$(OPEN) report.html

lcov: test
	lcov -t “test” -o lcov_test.info -c -d .
	genhtml -o report lcov_test.info
	$(OPEN) ./report/index.html


clean:
	rm -f *.out *.gcda *.gcno *.o *.a report.* s21_test test_with_gcov* *.gcov s21_test.c *.info
	rm -rf ./lcov_report report/ log.txt
