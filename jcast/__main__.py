# -*- coding: utf-8 -*-

"""jcast.__main__: executed when jcast directory is called as script."""

import sys


if "--gui" in sys.argv[1:]:
    from jcast.gui import launch_gui

    sys.exit(launch_gui())
else:
    from jcast.main import main

    main()
