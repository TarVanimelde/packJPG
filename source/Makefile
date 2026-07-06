# packXXX Makefile, based on UniMake: Universal Makefile
# Created by Matthias Stirner, 01/2016

TARGET   = packJPG
CC       = gcc
CPP      = g++
RC       = windres -O coff
CPPFLAGS = -I. -O3 -Wall -pedantic -funroll-loops -ffast-math -fsched-spec-load -fomit-frame-pointer -std=c++14
LDFLAGS  = -static -static-libgcc -static-libstdc++
CSRC     = $(wildcard *.c)
CPPSRC   = $(wildcard *.cpp)
DEPS     = $(wildcard *.h) Makefile
OBJ      = $(patsubst %.c,%.o,$(CSRC)) $(patsubst %.cpp,%.o,$(CPPSRC))

# conditional stuff
ifeq ($(OS),Windows_NT)
LDFLAGS  += -lpthread -L libwinpthread-1.dll
RES       = icons.res
UPX      := -upx --best --lzma $(TARGET).exe
else
CPPFLAGS += -DUNIX
RC        = 
RES       =
UPX       =
endif

%.o: %.cpp $(DEPS)
	$(CPP) -c -o $@ $< $(CPPFLAGS)
	
%.res: %.rc
	@-$(RC) $< $@

$(TARGET): $(OBJ) $(RES)
	$(CPP) -o $@ $^ -s $(LDFLAGS)
	$(UPX)

.PHONY: all dev lib dll test test-cli test-asan test-ubsan fixtures clean-test

all: $(TARGET)

dev: CPPFLAGS += -DDEV_BUILD
dev: $(TARGET)

lib: CPPFLAGS += -DBUILD_LIB
lib: $(OBJ)
	ar r $(TARGET)lib.a $(OBJ)
	ranlib $(TARGET)lib.a
    
dll: CPPFLAGS += -DBUILD_DLL
dll: LDFLAGS  += -Wl,--out-implib,libpackJPG.a -fvisibility=hidden
dll: $(OBJ)
	$(CPP) -shared -o $(TARGET).dll $^ $(LDFLAGS)

# -----------------------------------------------------------------------------
# Test suite (see tests/README.md)
#   make test       build + run the C++ unit/round-trip suite AND the CLI golden
#                   test against the real packJPG binary.
#   make fixtures   regenerate the JPEG fixture corpus (needs python3+Pillow and
#                   ImageMagick `convert`); the corpus is committed, so this is
#                   only needed when changing fixtures.
# -----------------------------------------------------------------------------
TESTDIR    = tests
TEST_SRC   = $(wildcard $(TESTDIR)/*.cpp)
TEST_ENG   = packjpg.cpp aricoder.cpp bitops.cpp
TEST_BIN   = $(TESTDIR)/run_tests
# BUILD_LIB compiles packjpg.cpp WITHOUT main() and exposes the pjglib_* API.
TESTFLAGS  = -I. -std=c++14 -O2 -g -DBUILD_LIB -DUNIX \
             -DFIXTURE_DIR=\"$(TESTDIR)/fixtures\"
TESTLD     =

$(TEST_BIN): $(TEST_SRC) $(TEST_ENG) $(wildcard $(TESTDIR)/*.h) $(DEPS)
	$(CPP) $(TESTFLAGS) $(TEST_SRC) $(TEST_ENG) -o $(TEST_BIN) $(TESTLD)

test: $(TEST_BIN) $(TARGET)
	@echo "==== C++ unit / round-trip suite ===="
	./$(TEST_BIN)
	@echo
	@echo "==== CLI golden test ===="
	@bash $(TESTDIR)/cli_golden.sh ./$(TARGET)

# CLI golden test only (broad corpus round-trip against the real binary).
test-cli: $(TARGET)
	@bash $(TESTDIR)/cli_golden.sh ./$(TARGET)

# AddressSanitizer + LeakSanitizer memory-safety regression guard. Builds packJPG
# with ASan (detect_leaks=1) and (a) runs a full compress+decompress round-trip
# over every VALID fixture and (b) feeds every INVALID fixture (malformed .pjg to
# the decoder, malformed .jpg to the encoder) to exercise both paths on hostile
# input. This guards the memory bugs fixed in this tree:
#   * the model_b heap-buffer-overflow (aricoder.cpp:629), root-caused to the
#     padbit char-signedness issue (packjpg.cpp:564);
#   * the new[]/free() alloc-dealloc mismatch (bitops.cpp get_c_bytes);
#   * the out-of-bounds heap write in pjg_unoptimize_header (packjpg.cpp), which
#     corrupted the allocator on malformed headers -> crash in free() (see
#     fixtures/invalid/pjg_header_oob.pjg);
#   * the qtables out-of-bounds index (jpg_setup_imginfo, see
#     fixtures/invalid/jpg_qtable_index_oob.jpg);
#   * leaks on error paths (read_jpeg `segment`, decode_jpeg `huffr`, unpack_pjg
#     `decoder`) reported by LeakSanitizer as issues #29 / #34;
#   * the errormessage[128] global-buffer-overflow from unbounded sprintf/strcpy
#     of long filenames / exception text (issue #30), guarded by feeding the CLI
#     an over-long filename below.
# For invalid inputs a non-zero exit is fine (rejection); only an actual ASan or
# LeakSanitizer report counts as a failure. It is EXPECTED TO PASS (exit 0); any
# such report means a regression. Requires a Linux ASan-capable g++/clang++
# (LeakSanitizer is not available on macOS). (UBSan is intentionally omitted --
# there is a separate, benign negative-shift UB documented in tests/README.md.)
ASANFLAGS = -std=c++14 -DUNIX -I. -O1 -g -fsanitize=address -fno-omit-frame-pointer
test-asan:
	$(CPP) $(ASANFLAGS) $(TEST_ENG) -o $(TESTDIR)/packJPG_asan $(TESTLD)
	@echo "Running packJPG under AddressSanitizer (compress+decompress, valid fixtures)..."
	@fail=0; for f in $(TESTDIR)/fixtures/valid/*.jpg; do \
	   d=$$(mktemp -d); cp "$$f" $$d/img.jpg; \
	   ASAN_OPTIONS=detect_leaks=1 ./$(TESTDIR)/packJPG_asan -p -np $$d/img.jpg >$$d/log 2>&1 || fail=1; \
	   rm -f $$d/img.jpg; \
	   ASAN_OPTIONS=detect_leaks=1 ./$(TESTDIR)/packJPG_asan -p -np $$d/img.pjg >>$$d/log 2>&1 || fail=1; \
	   if grep -qE "AddressSanitizer|LeakSanitizer" $$d/log; then echo "  ASAN-ERROR $$(basename $$f)"; fail=1; \
	   else echo "  clean      $$(basename $$f)"; fi; rm -rf $$d; \
	 done; \
	 echo "Feeding malformed fixtures (.pjg decompress + .jpg compress) to AddressSanitizer (rejection is fine, ASan report is not)..."; \
	 for f in $(TESTDIR)/fixtures/invalid/*.pjg $(TESTDIR)/fixtures/invalid/*.jpg; do \
	   d=$$(mktemp -d); ext=$${f##*.}; cp "$$f" $$d/bad.$$ext; \
	   ASAN_OPTIONS=detect_leaks=1 ./$(TESTDIR)/packJPG_asan -np $$d/bad.$$ext >$$d/log 2>&1 || true; \
	   if grep -qE "AddressSanitizer|LeakSanitizer" $$d/log; then echo "  ASAN-ERROR $$(basename $$f)"; fail=1; \
	   else echo "  clean      $$(basename $$f)"; fi; rm -rf $$d; \
	 done; \
	 echo "Feeding an over-long filename to the CLI error path under ASan (issue #30: errormessage[128] overflow)..."; \
	 d=$$(mktemp -d); longname="nonexistent_$$(head -c 200 /dev/zero | tr '\0' A)"; \
	 ASAN_OPTIONS=detect_leaks=1 ./$(TESTDIR)/packJPG_asan -np "$$d/$$longname" >$$d/log 2>&1 || true; \
	 if grep -qE "AddressSanitizer|LeakSanitizer" $$d/log; then echo "  ASAN-ERROR long-filename"; fail=1; \
	 else echo "  clean      long-filename (#30 guard)"; fi; rm -rf $$d; \
	 if [ $$fail -eq 0 ]; then echo "test-asan: PASS (no memory errors)"; else echo "test-asan: FAIL"; exit 1; fi

# UndefinedBehaviorSanitizer regression guard. Builds packJPG with UBSan and
# round-trips every valid fixture (compress + decompress). This guards the two
# UB classes fixed in this tree, both of which fired on ordinary valid images:
#   * negative shift exponent in the DEVLI macro (1 << (s-1) with s == 0);
#   * signed left shift of negative DCT coefficients (`x << cs_sal`).
# A UBSan "runtime error:" line means a regression. EXPECTED TO PASS (exit 0).
# Requires a UBSan-capable g++/clang++.
UBSANFLAGS = -std=c++14 -DUNIX -I. -O1 -g -fsanitize=undefined -fno-omit-frame-pointer
test-ubsan:
	$(CPP) $(UBSANFLAGS) $(TEST_ENG) -o $(TESTDIR)/packJPG_ubsan $(TESTLD)
	@echo "Round-tripping valid fixtures under UndefinedBehaviorSanitizer..."
	@fail=0; for f in $(TESTDIR)/fixtures/valid/*.jpg; do \
	   d=$$(mktemp -d); cp "$$f" $$d/img.jpg; \
	   UBSAN_OPTIONS=halt_on_error=0 ./$(TESTDIR)/packJPG_ubsan -p -np $$d/img.jpg >$$d/log 2>&1 || true; \
	   rm -f $$d/img.jpg; \
	   UBSAN_OPTIONS=halt_on_error=0 ./$(TESTDIR)/packJPG_ubsan -p -np $$d/img.pjg >>$$d/log 2>&1 || true; \
	   if grep -q "runtime error:" $$d/log; then echo "  UBSAN-ERROR $$(basename $$f)"; \
	     grep "runtime error:" $$d/log | head -2 | sed 's/^/    /'; fail=1; \
	   else echo "  clean      $$(basename $$f)"; fi; rm -rf $$d; \
	 done; \
	 if [ $$fail -eq 0 ]; then echo "test-ubsan: PASS (no undefined behavior)"; else echo "test-ubsan: FAIL"; exit 1; fi

fixtures:
	python3 $(TESTDIR)/gen_fixtures.py

clean-test:
	@-rm -f $(TEST_BIN)

clean: clean-test
	@echo clean...
	@-rm *.o *.a $(TARGET) $(TARGET).exe $(TARGET).dll
