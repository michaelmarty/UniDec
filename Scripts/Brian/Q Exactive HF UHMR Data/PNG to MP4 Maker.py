"""Combine the PNG images in C:\\Python\\Animations into an MP4 video."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterator

try:
    import imageio_ffmpeg
    from PIL import Image, ImageOps
except ImportError as exc:
    raise SystemExit(
        "Missing video dependencies. Install them with:\n"
        "    pip install Pillow imageio-ffmpeg"
    ) from exc


DEFAULT_INPUT_DIR = Path(r"C:\Python\Animations")
DEFAULT_OUTPUT_NAME = "animation.mp4"


def natural_sort_key(path: Path) -> tuple[object, ...]:
    """Return a key that sorts embedded numbers by value instead of as text."""
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", path.name)
    )


def find_png_files(folder: Path) -> list[Path]:
    """Find all PNG files directly inside *folder* in natural filename order."""
    if not folder.is_dir():
        raise FileNotFoundError(f"Input folder does not exist: {folder}")

    png_files = sorted(
        (
            path
            for path in folder.iterdir()
            if path.is_file() and path.suffix.lower() == ".png"
        ),
        key=natural_sort_key,
    )
    if not png_files:
        raise FileNotFoundError(f"No PNG files found in: {folder}")
    return png_files


def even_size(size: tuple[int, int]) -> tuple[int, int]:
    """Round dimensions up to even values required by YUV 4:2:0 video."""
    width, height = size
    return width + width % 2, height + height % 2


def frame_bytes(
    image_path: Path,
    reference_size: tuple[int, int],
    video_size: tuple[int, int],
) -> bytes:
    """Load one image and fit it onto an RGB frame without distortion."""
    with Image.open(image_path) as source:
        image = ImageOps.exif_transpose(source).convert("RGB")

        if image.size != reference_size:
            image = ImageOps.contain(
                image,
                video_size,
                method=Image.Resampling.LANCZOS,
            )

        if image.size != video_size:
            canvas = Image.new("RGB", video_size, "black")
            left = (video_size[0] - image.width) // 2
            top = (video_size[1] - image.height) // 2
            canvas.paste(image, (left, top))
            image = canvas

        return image.tobytes()


def progress_frames(
    png_files: list[Path],
    reference_size: tuple[int, int],
    video_size: tuple[int, int],
) -> Iterator[bytes]:
    """Yield encoded frame data while reporting progress."""
    total = len(png_files)
    report_every = max(1, total // 20)

    for index, image_path in enumerate(png_files, start=1):
        if index == 1 or index == total or index % report_every == 0:
            print(f"Processing frame {index}/{total}: {image_path.name}")
        yield frame_bytes(image_path, reference_size, video_size)


def create_video(input_dir: Path, output_file: Path, fps: float, quality: int) -> None:
    """Encode every PNG in *input_dir* as an H.264 MP4 video."""
    if fps <= 0:
        raise ValueError("FPS must be greater than zero.")
    if not 0 <= quality <= 51:
        raise ValueError("Quality must be between 0 (best) and 51 (smallest).")

    png_files = find_png_files(input_dir)
    output_file = output_file.resolve()
    output_file.parent.mkdir(parents=True, exist_ok=True)

    with Image.open(png_files[0]) as first_image:
        reference_size = first_image.size
    video_size = even_size(reference_size)

    print(f"Found {len(png_files)} PNG files in {input_dir.resolve()}")
    print(f"Creating {output_file} at {fps:g} FPS ({video_size[0]}x{video_size[1]})")

    writer = imageio_ffmpeg.write_frames(
        str(output_file),
        video_size,
        fps=fps,
        codec="libx264",
        pix_fmt_in="rgb24",
        pix_fmt_out="yuv420p",
        quality=None,
        macro_block_size=2,
        output_params=[
            "-crf",
            str(quality),
            "-preset",
            "medium",
            "-movflags",
            "+faststart",
        ],
    )

    writer.send(None)
    try:
        for frame in progress_frames(png_files, reference_size, video_size):
            writer.send(frame)
    finally:
        writer.close()

    print(f"Finished: {output_file}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine naturally sorted PNG files into an H.264 MP4 video."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"folder containing PNG files (default: {DEFAULT_INPUT_DIR})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=f"output MP4 path (default: <input-dir>\\{DEFAULT_OUTPUT_NAME})",
    )
    parser.add_argument(
        "--fps",
        type=float,
        default=60.0,
        help="video frame rate (default: 60)",
    )
    parser.add_argument(
        "--quality",
        type=int,
        default=18,
        help="H.264 CRF quality, 0-51; lower is better (default: 18)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_file = args.output or args.input_dir / DEFAULT_OUTPUT_NAME
    create_video(args.input_dir, output_file, args.fps, args.quality)


if __name__ == "__main__":
    main()
