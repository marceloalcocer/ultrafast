#!/usr/bin/env python

from distutils.core import setup

setup(
	name="ultrafast",
	version="0.1",
	description="Ultrafast optics package",
	author="Marcelo J P Alcocer",
	author_email="marcelo.j.p.alcocer@gmail.com",
	url="https://github.com/marceloalcocer/ultrafast",
	packages=["ultrafast"],
	requires=["pyyaml", "scipy"],
	provides=["ultrafast"]
)
