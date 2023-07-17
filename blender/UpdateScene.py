import bpy

def update_scene(scene):
    frame = scene.frame_current

    dens_coll = None
    for collection in bpy.data.collections:
        if "densities" in collection.name:
            dens_coll = collection

    for obj in dens_coll.objects:
        obj.hide_set(True)
        obj.hide_render=True
    
        if obj.name.startswith("dens_%03d" % frame):
            obj.hide_set(False)
            obj.hide_render=False

if __name__=="__main__":
    scene = bpy.context.scene
    
    update_scene(scene)
    bpy.app.handlers.frame_change_pre.append(update_scene)