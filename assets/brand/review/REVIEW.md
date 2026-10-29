# Visual review — 2026-10-29

**English** | [简体中文](REVIEW.zh-CN.md)

Followed `/work/test/MapPloid/codex读图提示词.txt`: read format, dimensions,
pixel count, color mode and file size before viewing; checked PNG integrity;
opened only one image per visual call. The metadata inventory records every PNG.

| Image | Dimensions | Pixels | Review method |
|---|---|---|---|
| Original generated concept | 1254 × 1254 | 1,572,516 | Direct, below 2 million pixels; no reduction (0%) |
| Final wordmark | 1160 × 400 | 464,000 | Direct, no reduction (0%) |
| Final opaque small icons | 32 × 32; 64 × 64 | 1,024; 4,096 | Each viewed at its exported size |
| Background comparison sheet | 1200 × 800 | 960,000 | One composite review image |

The comparison sheet uses 300 × 300 tiles derived from 1024 × 1024 transparent,
black and white icon exports: 91.42% fewer pixels per tile. Its white wordmark
tile is 540 × 186, derived from 1160 × 400: 78.35% fewer pixels. The composite is
not a uniformly downsampled copy of the original concept. All tile resampling
uses Lanczos, preserves proportions to the nearest pixel and leaves sources
unchanged. PNG encoding is lossless; resolution reduction is not strictly
lossless. The sheet is for layout and background review, not a replacement for
the original or final artwork. Its dimensions and integrity were verified
before visual inspection.

Completed: full concept and wordmark inspection; 32/64 px silhouette checks;
light/dark background checks for color, black and white versions; wordmark
spelling, spacing, clipping, alignment, color and clear-space inspection.
No local crop was needed; no separate crop was inspected. Other resolutions
were checked for integrity and metadata, not individually visually inspected.
No live GitHub browser rendering was tested.

Findings: the initial raster had colored edge speckles and a broad gold segment
(minor visual defects). Native SVG curves remove these artifacts and use the
specified colors. No clipping, overlap or severe visual defect was found in
the reviewed final assets. At 32 px, the gold accent is subdued (minor size
limitation); the M and baseline remain distinguishable. White marks require
a dark background, black marks a light background.

The original generated image and its local archived copy are byte-identical, verified
by SHA-256. The original PNG remains in the local archive and is not distributed here.
Original dimensions, files and pixels were not changed. Final SVGs
are separate artwork. See `image_metadata.json` for dimensions, bytes, color
modes and hashes, and `mapploid-backgrounds_visual_review_1200px.png` for the sheet.
