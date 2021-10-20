import { Dataset, DatasetType, IVector, Loader, useCancellablePromise } from "projection-space-explorer";
export declare class SDFLoader implements Loader {
    vectors: IVector[];
    datasetType: DatasetType;
    loading_area: string;
    resolvePath(entry: any, finished: (dataset: Dataset) => void, cancellablePromise?: ReturnType<typeof useCancellablePromise>["cancellablePromise"], modifiers?: string, abort_controller?: any): void;
    resolveContent(file: any, finished: (dataset: Dataset) => void, cancellablePromise?: ReturnType<typeof useCancellablePromise>["cancellablePromise"], modifiers?: string, controller?: any): void;
    loadCSV(finished: (dataset: Dataset) => void, entry: any, cancellablePromise?: ReturnType<typeof useCancellablePromise>["cancellablePromise"], modifiers?: string, controller?: any): void;
}
//# sourceMappingURL=SDFLoader.d.ts.map