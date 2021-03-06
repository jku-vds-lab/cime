import { Vect } from "../../Utility/Data/Vect";

/**
 * Helper class for cases where cross interaction could occur, one example would be that in the storytelling
 * the other lines should be grayed out, but a selection should overwrite this. This is solved by having 'layers'
 * of boolean values, each for one feature (1 layer for selection, 1 layer for storytelling) and solving the
 * actual value by the layer hierarchy.
 */
export class LayeringSystem {
    layers: Layer[] = []
    layerSize: number
    layerDictionary: { [index: number]: number } = {}
    
    constructor(layerSize: number) {
        this.layerSize = layerSize
    }

    clearLayer(layer: number, value: boolean) {
        let e = this.layers[this.layerDictionary[layer]]
        for (let i = 0; i < this.layerSize; i++) {
            e.setValue(i, value)
        }
    }

    setLayerActive(layer: number, active: boolean) {
        this.layers[this.layerDictionary[layer]].active = active
    }

    registerLayer(layer: number, active: boolean) {
        let i = this.layers.push(new Layer(this.layerSize, active)) - 1
        this.layerDictionary[layer] = i
    }

    setValue(index: number, layer: number, value: boolean) {
        let layerIndex = this.layerDictionary[layer]
        this.layers[layerIndex].setValue(index, value)
    }

    getValue(index: number) {
        let activeLayerIndex = this.layerDictionary[Object.keys(this.layerDictionary).sort((a, b) => Number.parseInt(b) - Number.parseInt(a)).find(key => this.layers[this.layerDictionary[key]].active)]

        //console.log(activeLayerIndex)

        if (activeLayerIndex == undefined || activeLayerIndex < 0) {
            return false
        }

        return this.layers[activeLayerIndex].getValue(index)
    }
}


class Layer {
    values: boolean[]
    active: boolean

    constructor(size: number, active: boolean) {
        this.values = new Array(size).fill(false)
        this.active = active
    }

    setValue(index: number, value: boolean) {
        this.values[index] = value
    }

    getValue(index: number) {
        return this.values[index]
    }
}

