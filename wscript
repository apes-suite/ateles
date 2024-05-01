#!/usr/bin/env python3
# encoding: utf-8
# 2023 Harald Klimach <harald.klimach@dlr.de>

APPNAME = 'ateles'

top = '.'
out = 'build'

def options(opt):
    '''Building options provided by Ateles.
       Remember, all options can be displayed with waf --help.'''
    opt.recurse('bin')
    opt.recurse('aotus')
    opt.recurse('tem')
    opt.recurse('atl')


def configure(conf):
    '''Project configuration'''
    import os
    conf.recurse('aotus', 'subconf')
    conf.recurse('bin', 'preconfigure')
    # Initialize the coco preprocessing tool
    conf.load('coco')
    conf.env['COCOSET'] = 'default.coco'
    if not conf.options.coco_reports:
      # Make coco silent, if not explicitly asked for reports:
      if conf.env.COCOFLAGS:
        conf.env.COCOFLAGS.insert(0, '-s')
        conf.env.COCOFLAGS.append('-ad')
      else:
        conf.env.COCOFLAGS = ['-s', '-ad']
    conf.recurse('tem')
    conf.recurse('polynomials')
    conf.recurse('atl')
    conf.recurse('bin', 'postconfigure')


def build(bld):
    '''Build the Ateles project'''
    from revision_module import fill_revision_string
    bld.recurse('bin')
    if not (bld.cmd == 'docu' and bld.env.fordonline):
        bld.recurse('aotus')
    fill_revision_string(bld, subdir='atl')
    bld(rule='cp ${SRC} ${TGT}', source=bld.env.COCOSET, target='coco.set')
    if not (bld.cmd == 'docu' and bld.env.fordonline):
        bld.recurse('tem')
    else:
        bld(rule='cp ${SRC} ${TGT}',
            source = bld.path.find_node(['tem', 'source', 'arrayMacros.inc']),
            target = bld.path.find_or_declare('arrayMacros.inc'))
        bld.add_group()
    bld.recurse('polynomials')
    bld.recurse('atl')

    if not (bld.cmd == 'docu'):
        bld.add_group()

        bld(
            rule = 'cp ${SRC} ${TGT[0].abspath()}',
            source = bld.path.find_or_declare('atl/ateles'),
            target = 'ateles'
        )
        bld(
            rule = 'cp ${SRC} ${TGT[0].abspath()}',
            source = bld.path.find_or_declare('atl/atl_harvesting'),
            target = 'atl_harvesting'
        )
        bld(
            rule = 'cp ${SRC} ${TGT[0].abspath()}',
            source = bld.path.find_or_declare('aotus/lua'),
            target = 'lua'
        )

#clean build directory and coco completely to create the build from scratch
def cleanall(ctx):
    from waflib import Options
    Options.commands = ['distclean'] + Options.commands
    ctx.exec_command('rm coco')
