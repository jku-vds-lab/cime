import { format } from "d3";
import { ScaleMappingFunction, } from "lineupjs";
export var DEFAULT_FORMATTER = format(".3n");
export function noNumberFilter() {
    return { min: -Infinity, max: Infinity, filterMissing: false };
}
export function isEqualNumberFilter(a, b, delta) {
    if (delta === void 0) { delta = 0.001; }
    return (similar(a.min, b.min, delta) &&
        similar(a.max, b.max, delta) &&
        a.filterMissing === b.filterMissing);
}
export function similar(a, b, delta) {
    if (delta === void 0) { delta = 0.5; }
    if (a === b) {
        return true;
    }
    return Math.abs(a - b) < delta;
}
export function isUnknown(v) {
    return v == null || v === undefined || isNaN(v);
}
export function isDummyNumberFilter(filter) {
    return (!filter.filterMissing && !isFinite(filter.min) && !isFinite(filter.max));
}
export function restoreMapping(desc, factory) {
    if (desc.map) {
        return factory.mappingFunction(desc.map);
    }
    return new ScaleMappingFunction(desc.domain || [0, 1], "linear", desc.range || [0, 1]);
}
export function restoreNumberFilter(v) {
    return {
        min: v.min != null && isFinite(v.min) ? v.min : -Infinity,
        max: v.max != null && isFinite(v.max) ? v.max : +Infinity,
        filterMissing: v.filterMissing,
    };
}
//# sourceMappingURL=helper_methods.js.map