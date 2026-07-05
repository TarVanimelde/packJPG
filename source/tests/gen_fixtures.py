#!/usr/bin/env python3
"""Generate the JPEG fixture corpus for packJPG round-trip tests.

packJPG ships no real .jpg samples (docs/sample_images.zip contains only PNG
coefficient dumps), so the round-trip suite needs its own corpus.

IMPORTANT encoder note (discovered empirically, see tests/README.md):
  * packJPG is picky about JPEG Huffman *pad bits* (packjpg.cpp:2726-2732,
    "inconsistent use of padbits"). Pillow's JPEG encoder pads with 0-bits and
    packJPG rejects those outright. libjpeg (via ImageMagick `convert`) pads per
    the JPEG spec, and packJPG accepts those files -- so we encode with libjpeg.
  * Even with libjpeg, many files still raise the padbit warning (it depends on
    where the Huffman bitstream happens to end), which packJPG treats as an
    error under DEFAULT options. Passing `-p` (proceed-on-warnings) makes packJPG
    use padbit=1 -- exactly what libjpeg wrote -- so the round-trip stays
    *bit-identical*. Whether a given file needs `-p` is content-dependent, so the
    in-process test probes acceptance at runtime and the CLI golden test simply
    uses `-p` across the whole corpus.

Corpus layout:
  fixtures/valid/    libjpeg JPEGs (grayscale/color/progressive/cmyk/odd/...);
                     all round-trip bit-identically with `-p`.
  fixtures/invalid/  malformed inputs for error-handling tests.

We generate raw pixels with Pillow (deterministic) but ENCODE to JPEG with
ImageMagick's `convert` (libjpeg). Run from source/tests/:
    python3 gen_fixtures.py
Deterministic and idempotent.
"""
import os
import subprocess
import sys

from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")
VALID = os.path.join(FIX, "valid")
INVALID = os.path.join(FIX, "invalid")


def content_png(path, w, h, mode="RGB"):
    """Deterministic image with gradients + shapes so DCT blocks are non-trivial."""
    img = Image.new(mode, (w, h))
    px = img.load()
    for y in range(h):
        for x in range(w):
            r = (x * 255) // max(1, w - 1)
            g = (y * 255) // max(1, h - 1)
            b = ((x + y) * 255) // max(1, w + h - 2)
            px[x, y] = (r + g + b) // 3 if mode == "L" else (r, g, b)
    d = ImageDraw.Draw(img)
    d.ellipse([w // 6, h // 6, w * 5 // 6, h * 5 // 6],
              outline=(0 if mode == "L" else (0, 0, 0)), width=2)
    d.rectangle([w // 3, h // 3, w * 2 // 3, h * 2 // 3],
                outline=(255 if mode == "L" else (255, 255, 255)), width=1)
    img.save(path)


def convert(src_png, dst_jpg, *opts):
    subprocess.run(["convert", src_png, "-strip", *opts, dst_jpg],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(f"    {os.path.relpath(dst_jpg, FIX):40s} {os.path.getsize(dst_jpg):>7d} bytes")


def main():
    for d in (VALID, INVALID):
        os.makedirs(d, exist_ok=True)
    tmp = os.path.join(FIX, "_src")
    os.makedirs(tmp, exist_ok=True)

    def png(name, w, h, mode="RGB"):
        p = os.path.join(tmp, name)
        content_png(p, w, h, mode)
        return p

    print("  valid/ (libjpeg JPEGs, round-trip bit-identical with -p):")
    convert(png("g40.png", 40, 40, "L"), os.path.join(VALID, "grayscale.jpg"),
            "-quality", "80")
    convert(png("g128.png", 128, 96, "L"), os.path.join(VALID, "grayscale_large.jpg"),
            "-quality", "88")
    c = png("c64.png", 64, 64)
    convert(c, os.path.join(VALID, "baseline_rgb.jpg"), "-quality", "85",
            "-sampling-factor", "2x2")
    convert(png("c48.png", 48, 48), os.path.join(VALID, "baseline_444.jpg"),
            "-quality", "95", "-sampling-factor", "1x1")
    convert(c, os.path.join(VALID, "baseline_q20.jpg"), "-quality", "20")
    convert(c, os.path.join(VALID, "progressive.jpg"), "-quality", "85",
            "-interlace", "JPEG")
    convert(png("c32.png", 32, 32), os.path.join(VALID, "cmyk.jpg"), "-quality", "85",
            "-colorspace", "CMYK")
    convert(png("odd.png", 17, 13), os.path.join(VALID, "odd_dimensions.jpg"),
            "-quality", "85")
    convert(png("c256.png", 256, 192), os.path.join(VALID, "large_256.jpg"),
            "-quality", "88")

    print("  invalid/ (malformed inputs):")

    def winvalid(name, data):
        p = os.path.join(INVALID, name)
        with open(p, "wb") as f:
            f.write(data)
        print(f"    invalid/{name:32s} {len(data):>7d} bytes")

    winvalid("empty.bin", b"")
    winvalid("not_a_jpeg.txt", b"this is plainly not an image file\n" * 4)
    winvalid("bad_after_soi.jpg", b"\xFF\xD8" + b"\x00\x11\x22\x33" * 8)
    with open(os.path.join(VALID, "large_256.jpg"), "rb") as f:
        full = f.read()
    winvalid("truncated.jpg", full[:200])
    winvalid("soi_eoi_only.jpg", b"\xFF\xD8\xFF\xD9")
    # Starts with the PJG magic ('JS', packjpg.cpp:711) but is otherwise garbage,
    # exercising the decompress path's rejection of a corrupt PJG.
    winvalid("fake_pjg.pjg", b"JS" + b"\x00\xFF" * 16)

    # tidy up intermediate PNGs
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)
    print("\n  Fixture generation complete.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:  # noqa
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
