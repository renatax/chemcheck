from django.shortcuts import render
# Create your views here.
from django.shortcuts import render_to_response
from models import EntryForm
from django.template import Template
from django.template import RequestContext
#import json
from django.contrib.auth.decorators import login_required

@login_required

def store(request):
    if request.method == 'POST':
            name = request.POST['name']
            reaction = request.POST['reaction']
            products = request.POST['products']
            conditions = request.POST['conditions']
            synthons = request.POST['synthons']
            doi = request.POST['doi']
            
            try:
                #now everything should be ok, so we append write the reaction on a file...
                #filename="reactions.txt"
                #f=open(filename, 'a')
#                 f=write(name)
#                 f.write("\t"+reaction)
#                 f.write("\t"+products)
#                 f.write("\t"+conditions)
#                 f.write("\t"+synthons)
#                 f.write("\t"+doi+"\n")
#                 f.close()
                reactionstring=name + "\t" + reaction +"\t" + products + "\t" + conditions + "\t" + str(synthons)+"\t"+doi

            except:
                pass
    return render_to_response('templates/results.html', {'reactionstring':reactionstring},context_instance=RequestContext(request))


@login_required
def enter(request):
    def errorHandle(error):
        return render_to_response('templates/default.html', {
                'error' : error,
                'form' : form,
        },context_instance=RequestContext(request))
        
    if request.method == 'POST': # If the form has been submitted...
        form = EntryForm(request.POST)
        reaction = request.POST['reaction']
        out = form.getOut()
        stuff=form.getStuff()
        
        fullreaction=form.getFullReactions()
        print fullreaction
        if form.is_valid():# and formisreallyok: # All validation rules pass
            name = request.POST['name']
            reaction = request.POST['reaction']
            products = request.POST['products']
            conditions = request.POST['conditions']
            synthons = request.POST['synthons']
            doi = request.POST['doi']
            #sadly /t and other tabulations get collapsed to one space in HTML. Therefore the fields are separated by 4*&nbsp;
            reactionstring=name + 4*"&nbsp;" + reaction +4*"&nbsp;" + products + 4*"&nbsp;" + conditions + 4*"&nbsp;" + str(synthons)+4*"&nbsp;"+doi
            
            return render_to_response('templates/results.html', {
                        'form': form,
                        'reactionstring': reactionstring
                    },context_instance=RequestContext(request))
        else:
            error = u'form is invalid'
            return errorHandle(error)
    else:
        form = EntryForm() 
        
        out = form.getOut()# An unbound form
        return render_to_response('templates/default.html', {
            'form': form,
        }, context_instance=RequestContext(request))
