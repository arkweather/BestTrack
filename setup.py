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
			'shapely~=1.15.16',
			'scipy~=0.18.0',
			'basemap~=1.0.7']
			
if __name__ == '__main__':
	
	pkg_description = "Best Track is approximately equivalent to the Warning Decision Support System - \
						Integrated Information’s (WDSS-II) w2besttrack algorithm with the potential for \
						additional features and greater flexibility."

	setup(name="besttrack",
          version="0.1.0",
          description="Object-based geotemporal path optimization",
          author="David Harrison",
          author_email="david.r.harrison-1@ou.edu",
          long_description=pkg_description,
          license="MIT",
          url="https://github.com/arkweather/BestTrack",
          packages=["besttrack"],
          scripts=["bin/best_track"],
          keywords=["verification", "tracking", "weather", "meteorology"],
          classifiers=classifiers,
          install_requires=requires)
