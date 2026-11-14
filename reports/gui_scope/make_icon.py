# -*- coding: utf-8 -*-

"""The GUI's icon (PR 21; reports/2026-11-13-gui-cli-scope.md §13): the website's fish logo, its transparent margins
trimmed, centred on a transparent square and scaled to 512 x 512 pixels. Needs Pillow (a dependency of matplotlib).

Usage (from the repository root):
    python reports/gui_scope/make_icon.py [--logo jcast_website/images/fish.png] [--out jcast/gui/icon.png]
"""

import argparse

from PIL import Image

SIZE = 512


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--logo', default='jcast_website/images/fish.png')
    parser.add_argument('--out', default='jcast/gui/icon.png')
    args = parser.parse_args()

    logo = Image.open(args.logo).convert('RGBA')
    logo = logo.crop(logo.getchannel('A').getbbox())
    side = max(logo.size)
    square = Image.new('RGBA', (side, side), (0, 0, 0, 0))
    square.paste(logo, ((side - logo.width) // 2, (side - logo.height) // 2))
    square.resize((SIZE, SIZE), Image.LANCZOS).save(args.out, optimize=True)
    print(f'{args.logo} {Image.open(args.logo).size}: trimmed to {logo.size}, saved as {args.out} ({SIZE} x {SIZE})')


if __name__ == '__main__':
    main()
