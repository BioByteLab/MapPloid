# MapPloid identity

**English** | [简体中文](README.zh-CN.md)

![MapPloid](mapploid-wordmark.svg)

Approved concept A: two rounded frequency peaks form an M above a short baseline.
The mark alludes to allele-frequency modelling; it does not depict a particular
ploidy or a measured distribution. Colors: teal `#087F8C`, gold `#D89B45`,
off-white `#F7FAF9`. Keep the name exactly `MapPloid`.

| Use | SVG master | PNG |
|---|---|---|
| Square project icon | `mapploid-icon.svg` | `mapploid-icon-512.png` |
| Transparent icon | `mapploid-icon-transparent.svg` | `mapploid-icon-transparent-512.png` |
| README wordmark | `mapploid-wordmark.svg` | `mapploid-wordmark.png` |
| Transparent wordmark | `mapploid-wordmark-transparent.svg` | `mapploid-wordmark-transparent.png` |
| Black on light backgrounds | `mapploid-icon-black.svg`, `mapploid-wordmark-black.svg` | Matching PNG files |
| White on dark backgrounds | `mapploid-icon-white.svg`, `mapploid-wordmark-white.svg` | Matching PNG files |

All four icons have 32, 64, 256, 512 and 1024 px PNG exports. All four wordmarks
have 1160 × 400 px PNG exports. SVGs contain editable paths, not embedded raster
images. Lettering is outlined from DejaVu Sans Bold, so viewers need no font.
At 32 px the gold accent is subtle; the silhouette carries recognition.
Preserve the aspect ratio and clear space. Prefer the white version on dark
backgrounds and the opaque square icon for avatars.

Both project READMEs use theme-aware SVG wordmarks. Assets are included in the GitHub repository.

The initial concept used the built-in imagegen tool, with the exact
[prompt](source/imagegen-prompt.txt). Its unmodified PNG is retained in the original local project archive; the public
repository distributes the finished SVG/PNG artwork. Final SVG artwork was redrawn as native
curves to remove raster fringes, standardize colors, shorten the accent and
center the mark. It is a vector interpretation, not a pixel-identical trace.

Regenerate final assets with `python assets/brand/build.py` from the version
directory (requires matplotlib and ImageMagick `convert`). The script overwrites
the named SVG/PNG exports, while preserving the original concept. It is an asset
build helper, not a runtime dependency of MapPloid.

See [visual review](review/REVIEW.md) and [image metadata](review/image_metadata.json).
