from setuptools import setup
setup(
    name = 'averagemaxz',
    packages = ['averagemaxz'],
    entry_points = {
        'console_scripts': [
            'maxZ = averagemaxz.__main__:main',
            'maxz = averagemaxz.__main__:main'
        ]
    })
