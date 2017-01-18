#!/usr/bin/python
# -*- coding: utf-8 -*-

from setuptools import setup

classifiers = ['Development Status :: 4 - Beta',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: MIT License',
               'Programming Language :: Python :: 2',
               'Programming Language :: Python :: 2.6',
               'Programming Language :: Python :: 2.7',
              ]
              
requires = ['matplotlib~=1.5.1',
			'numpy~=1.11.1',
			'beautifulsoup4~=4.4.1',
			'shapely~=1.5.16',
			'scipy~=0.18.0',
			'basemap~=1.0.7']
			
if __name__ == '__main__':
	
	pkg_description = "Best Track is approximately equivalent to the Warning Decision Support System - \
Integrated Information (WDSS-II) w2besttrack algorithm with the potential for \
additional features and greater flexibility."
	
	#pkg_description = "test"					

	setup(name="besttrack",
          version="0.2.0",
          description="Object-based geotemporal path optimization",
          author="David Harrison",
          author_email="david.r.harrison-1@ou.edu",
          long_description=pkg_description,
          license="MIT",
          url="https://github.com/arkweather/BestTrack",
          packages=["besttrack"],
          include_package_data = True,
          scripts=["bin/besttrack", "bin/besttrack.config"],
          keywords=["verification", "tracking", "weather", "meteorology"],
          classifiers=classifiers,
          install_requires=requires)
