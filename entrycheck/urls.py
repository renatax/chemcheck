'''
Created on Mar 22, 2013

@author: ernia
'''
from django.conf.urls import patterns, include, url

from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'damnyoudjango.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),
    
    url(r'^enter/$', 'entrycheck.views.enter'),
    url(r'^results/$', 'entrycheck.views.store'),
    
)
