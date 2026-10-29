"""Build the approved peak-M vector identity and PNG exports.

Requires matplotlib (outlined wordmark) and ImageMagick convert (SVG rendering).
The AI concept remains untouched in source/; final artwork is native SVG.
"""
from pathlib import Path
import subprocess
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath

ROOT = Path(__file__).resolve().parent
TEAL, GOLD, PAPER = '#087F8C', '#D89B45', '#F7FAF9'
CURVE = ('M72 328 C136 312 122 148 174 148 '
         'C214 148 212 304 256 304 C300 304 298 148 338 148 '
         'C390 148 376 312 440 328')
# A subcurve of the right peak, obtained by splitting its cubic at t=0.36.
ACCENT = 'M338 148 C356.72 148 366.8864 169.2544 374.953472 197.206528'


def mark(color=TEAL, accent=GOLD):
    return (f'<g fill="none" stroke="{color}" stroke-width="36" '
            f'stroke-linecap="round" stroke-linejoin="round">'
            f'<path d="{CURVE}"/>'
            f'<path d="M120 380 H392"/>'
            + (f'<path d="{ACCENT}" stroke="{accent}" '
               'stroke-linecap="butt"/>' if accent else '') + '</g>')


def wordmark(color):
    tp = TextPath((0, 0), 'MapPloid', size=142,
                  prop=FontProperties(family='DejaVu Sans', weight='bold'))
    commands = []
    for vertices, code in tp.iter_segments():
        if code == 79:
            commands.append('Z')
        else:
            commands.append({1: 'M', 2: 'L', 3: 'Q', 4: 'C'}[code] +
                            ' '.join(f'{v:.5f}' for v in vertices))
    return (f'<path fill="{color}" transform="translate(380 236) scale(1 -1)" '
            f'd="{" ".join(commands)}"/>')


def svg(name, width, height, body):
    path = ROOT / f'{name}.svg'
    path.write_text(f'<svg xmlns="http://www.w3.org/2000/svg" '
                    f'width="{width}" height="{height}" '
                    f'viewBox="0 0 {width} {height}" role="img" '
                    f'aria-label="MapPloid">\n<title>MapPloid</title>\n'
                    f'{body}\n</svg>\n')
    return path


def render(source, width, destination):
    subprocess.run(['convert', '-background', 'none', str(source),
                    '-resize', str(width), str(destination)], check=True)


def main():
    for name, color, accent, bg in [
        ('mapploid-icon', TEAL, GOLD, PAPER),
        ('mapploid-icon-transparent', TEAL, GOLD, None),
        ('mapploid-icon-black', '#000000', None, None),
        ('mapploid-icon-white', '#FFFFFF', None, None),
    ]:
        body = (f'<rect width="512" height="512" fill="{bg}"/>' if bg else '')
        source = svg(name, 512, 512, body + mark(color, accent))
        for size in (32, 64, 256, 512, 1024):
            render(source, size, ROOT / f'{name}-{size}.png')
    for name, color, accent, bg in [
        ('mapploid-wordmark', TEAL, GOLD, PAPER),
        ('mapploid-wordmark-transparent', TEAL, GOLD, None),
        ('mapploid-wordmark-black', '#000000', None, None),
        ('mapploid-wordmark-white', '#FFFFFF', None, None),
    ]:
        body = (f'<rect width="1160" height="400" fill="{bg}"/>' if bg else '')
        body += f'<g transform="translate(16 30) scale(.66)">{mark(color, accent)}</g>'
        body += wordmark(color)
        source = svg(name, 1160, 400, body)
        render(source, 1160, ROOT / f'{name}.png')


if __name__ == '__main__':
    main()
