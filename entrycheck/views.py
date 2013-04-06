
# Create your views here.
from django.shortcuts import render_to_response
from models import EntryForm
from django.template import RequestContext
#import json


# def marvin(request):
#     return ('file:///Users/ernia/Documents/workspace/damnyoudjango/marvin/')
def enter(request):
    def errorHandle(error):
        #form = EntryForm()
        return render_to_response('templates/default.html', {
                'error' : error,
                'form' : form,
                #'reactionsmarts' : form.js_data,
        },context_instance=RequestContext(request))
        
    if request.method == 'POST': # If the form has been submitted...
        form = EntryForm(request.POST)
        reaction = request.POST['reaction']
        #js_data = json.dumps({'reactionsmarts': reaction})
#         print js_data
        out = form.getOut()
        stuff=form.getStuff()
        #try:
        fullreaction=form.getFullReactions()
        print fullreaction
#         if stuff.has_key('reactants'): 
#             fullreaction='.'.join(stuff['reactants']) + ">>" + '.'.join(stuff['products'])
#         else:
#             fullreaction=""
            
        if form.is_valid():# and formisreallyok: # All validation rules pass
            name = request.POST['name']
            reaction = request.POST['reaction']
            products = request.POST['products']
            conditions = request.POST['conditions']
            synthons = request.POST['synthons']
            doi = request.POST['doi']
            
            
            return render_to_response('templates/results.html', {
                        'form': form,
                    },context_instance=RequestContext(request))
        else:
            error = u'form is invalid'
            return errorHandle(error)
    else:
        form = EntryForm() 
        #js_data=json.dumps({"reactionsmarts":"CCOCC"})
        #print js_data
        out = form.getOut()# An unbound form
        return render_to_response('templates/default.html', {
            'form': form,
            #'reactionsmarts' : form.js_data,
        }, context_instance=RequestContext(request))
